#ifndef INDEX_BUILD_HPP
#define INDEX_BUILD_HPP

#include "../../bundled/Minimizers/lib/BreiZHMinimizer.hpp"
#include "../constants.hpp"
#include "../colorsets.hpp"


namespace kaminari {
namespace minimizer {

uint64_t index::get_file_size(const std::string& filen) const {
    struct stat file_status;
    if (stat(filen.c_str(), &file_status) < 0) {
        return -1;
    }
    return file_status.st_size;
}

void
index::build(build::options_t& build_parameters)
{
    //STEP 1 : PARSE FILES =====================================================
    if (build_parameters.verbose > 0) std::cerr << "Step 1 & 2: parsing " << nb_docs << " files and merging results\n";

    auto start_time = std::chrono::high_resolution_clock::now();

    std::string Bzhminmer_tmp_file = build_parameters.output_dirname + "/tmp/result";

    generate_minimizers(
        build_parameters.input_filenames, 
        Bzhminmer_tmp_file,
        build_parameters.output_dirname + "/tmp", 
        build_parameters.nthreads, 
        build_parameters.max_ram_MB, 
        build_parameters.k, 
        build_parameters.m, 
        8,
        "crumsort",
        build_parameters.verbose,
        false, false, false // -> dont skip minmer step, dont keep tmp files
    );

    vector<std::string>().swap(build_parameters.input_filenames); //free memory

    std::cerr << "Time for reading minimizers + sort them: " <<
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start_time).count() <<
        " milliseconds\n";

    std::string Bzhminmer_file = Bzhminmer_tmp_file + "." + std::to_string(nb_docs) + "c";

    // creating input reading buffer
    uint64_t elem_words = 1 + ((nb_docs + 63)/64); // number of uint64_t for 1 minmer + colors
    uint64_t buffer_bytes = 2 * constants::MB;
    uint64_t input_buffer_elems = buffer_bytes / (elem_words * sizeof(uint64_t)); // number of elems for a 2MB buffer 
    std::vector<uint64_t> input_buffer(input_buffer_elems * elem_words);

    FILE* finput = fopen(Bzhminmer_file.c_str(), "rb");
    if (!finput) throw std::runtime_error("Cannot open parsing-result file");

    ankerl::unordered_dense::set<uint64_t> unique_minmers;

    while (true) {
        size_t got = fread(input_buffer.data(), elem_words * sizeof(uint64_t), input_buffer_elems, finput);
        if (got == 0) break;
        for (size_t i = 0; i < got; ++i) {
            uint64_t minimizer = input_buffer[i * elem_words];
            unique_minmers.insert(minimizer);
        }
    }

    rewind(finput);

    uint64_t nb_elems = unique_minmers.size();

    //STEP 3 : BUILDING MPHF ===================================================
    {
        if (build_parameters.verbose > 0) std::cerr << "Step 3: building the MPHF for " << unique_minmers.size() << " minimizers\n";
        start_time = std::chrono::high_resolution_clock::now();
        int backup, redirect;
        fflush(stdout);
        backup = dup(1);
        redirect = open("/dev/null", O_WRONLY);
        dup2(redirect, 1);
        close(redirect);

        //MPHF via PTHash, n minmers, hf : minmer(uint64) -> [0, n-1]
        hf.build_in_internal_memory(unique_minmers.begin(), unique_minmers.size(), get_pthash_options(build_parameters));

        fflush(stdout);
        dup2(backup, 1);
        close(backup);
        assert(hf.num_keys() == unique_minmers.size());

        std::cerr << "Time for MPHF build: " << 
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - start_time).count() << 
            " milliseconds\n";

        std::cerr << "size of MPHF : " << hf.num_bits()/8 << " Bytes\n";

        ankerl::unordered_dense::set<uint64_t>().swap(unique_minmers); //free memory
    }

    // common variables for CM & CS
    uint64_t ram_budget = static_cast<uint64_t>(build_parameters.max_ram_MB * 0.75); //75% of allocated RAM for CM & CS

    std::string minmer_to_cid_tmp_file = build_parameters.output_dirname + "/tmp/minmer_cid.bin";
    uint64_t output_buffer_elems = buffer_bytes / (2 * sizeof(uint64_t)); // number of (minmer-cid) for a 2MB buffer 

    uint64_t cid = -1; // starts at -1 so first unique set becomes 0


    //STEP 4.1 : COLORSETS DEDUP + STORING ====================================
    {
        if (build_parameters.verbose > 0) std::cerr << "Step 4.1: colors deduplication and storing: ColorSets\n";
        start_time = std::chrono::high_resolution_clock::now();
        
        
        FILE* fout = fopen(minmer_to_cid_tmp_file.c_str(), "wb+");
        if (!fout) throw std::runtime_error("Cannot open tmp minmer-cid file");

        std::vector<uint64_t> output_buffer(output_buffer_elems * 2);

        uint64_t nb_words_colors = elem_words - 1;
        uint64_t* last_color_words = new uint64_t[nb_words_colors]();

        colorsets::builder cs_builder(
            build_parameters.output_dirname + "/colorsets", 
            nb_docs, 
            ram_budget * constants::MB
        );
        uint64_t output_pos = 0; // how many uint64_t currently in output buffer
        uint64_t processed = 0;  // number of elements processed so far
        

        while (processed < nb_elems){
            // Read up to input_buf_elems elements
            size_t to_read = std::min(input_buffer_elems, nb_elems - processed);
            size_t got = fread(input_buffer.data(), elem_words * sizeof(uint64_t), to_read, finput);
            if (got != to_read) {
                throw std::runtime_error("Error reading input file");
            }

            for (size_t j = 0; j < to_read; j++) {
                const uint64_t* elem = &input_buffer[j * elem_words];
                uint64_t key = elem[0];
                const uint64_t* color_words = elem + 1;

                if (processed == 0 || std::memcmp(color_words, last_color_words, nb_words_colors * sizeof(uint64_t)) != 0) {
                    std::vector<uint32_t> colors;
                    for (size_t w = 0; w < nb_words_colors; w++) {
                        uint64_t word = color_words[w];
                        for (size_t bit = 0; bit < 64; bit++) {
                            if (word & (1ULL << bit)) {
                                colors.push_back(static_cast<uint32_t>(w * 64 + bit));
                            }
                        }
                    }
                    cs_builder.add_color_set(colors.data(), colors.size());
                    cid++;
                    std::memcpy(last_color_words, color_words, nb_words_colors * sizeof(uint64_t));
                }

                // Output element: [key, cid+check_bits]
                uint64_t minmer_last_bits = key & ((1UL << build_parameters.b)-1);
                output_buffer[output_pos++] = hf(key); //hash of minmer
                output_buffer[output_pos++] = (cid << build_parameters.b) | minmer_last_bits;

                if (output_pos == output_buffer.size()) {
                    fwrite(output_buffer.data(), sizeof(uint64_t), output_pos, fout);
                    output_pos = 0;
                }
                processed++;
            }
        }

        // Flush remaining output
        if (output_pos > 0) {
            fwrite(output_buffer.data(), sizeof(uint64_t), output_pos, fout);
        }

        cs_builder.build();

        delete[] last_color_words;
        fclose(finput);
        fflush(fout);
        fclose(fout);
        
        std::cerr << "Time for ColorSets: " << 
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start_time).count() << 
        " milliseconds\n";
    }

    //STEP 5 : COLORMAPPER  =============================================
    {
        if (build_parameters.verbose > 0) std::cerr << "Step 4.2: mapping hash -> minimizer: ColorMapper\n";
        start_time = std::chrono::high_resolution_clock::now();

        FILE* ffinal = fopen(minmer_to_cid_tmp_file.c_str(), "rb");
        if (!ffinal) throw std::runtime_error("Cannot open tmp minmer-cid file");

        colormapper::builder cm_builder(
            build_parameters.output_dirname + "/colormapper", 
            nb_elems, 
            ceil(log2(cid)) + build_parameters.b, //bits for cid + bits for checking alien minmers 
            ram_budget * constants::MB
        );
        uint64_t high = cm_builder.init_chunk();
        uint64_t low = 0;
        
        std::vector<uint64_t> final_buffer(output_buffer_elems * 2);

        // Each pass goes through the file completely
        while (high < nb_elems) {
            size_t processed = 0;
            while (processed < nb_elems) {
                size_t to_read = std::min(output_buffer_elems, nb_elems - processed);
                size_t got = fread(final_buffer.data(), 2 * sizeof(uint64_t), to_read, ffinal);
                if (got != to_read) {
                    throw std::runtime_error("Error reading tmp minmer-cid file");
                }

                for (size_t j = 0; j < got; j++) {
                    uint64_t key = final_buffer[j * 2];
                    uint64_t val = final_buffer[j * 2 + 1];
                    if (key >= low && key < high) {
                        cm_builder.set(key, val);
                    }
                }

                processed += got;
            }

            rewind(ffinal);
            low = high;
            high += cm_builder.flush();
        }

        // Final pass for the last chunk
        size_t processed = 0;
        while (processed < nb_elems) {
            size_t to_read = std::min(output_buffer_elems, nb_elems - processed);
            size_t got = fread(final_buffer.data(), 2 * sizeof(uint64_t), to_read, ffinal);
            if (got != to_read) {
                throw std::runtime_error("Error reading tmp minmer-cid file");
            }

            for (size_t j = 0; j < got; j++) {
                uint64_t key = final_buffer[j * 2];
                uint64_t val = final_buffer[j * 2 + 1];
                if (key >= low && key < high) {
                    cm_builder.set(key, val);
                }
            }

            processed += got;
        }

        cm_builder.build();
        
        fclose(ffinal);

        std::cerr << "Time for ColorMapper: " << 
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start_time).count() << 
        " milliseconds\n";
    }

    std::remove(Bzhminmer_file.c_str());
    std::remove(minmer_to_cid_tmp_file.c_str());
}

} // namespace minimizer
} // namespace kaminari


#endif // INDEX_BUILD_HPP