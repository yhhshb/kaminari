#ifndef INDEX_BUILD_HPP
#define INDEX_BUILD_HPP

#include "../../bundled/Minimizers/lib/BreiZHMinimizer.hpp"
#include "../constants.hpp"
#include "../colorsets.hpp"
#include "../utils.hpp"


namespace kaminari {
namespace minimizer {

void
index::build(build::options_t& build_parameters)
{
    //STEP 1 : PARSE FILES =====================================================
    if (build_parameters.verbose >= 1) std::cout << "[I] Step 1 & 2: parsing " << nb_docs << " files and merging results\n";

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

    if (build_parameters.verbose >= 1){
        std::cout << "[I] Step 1 & 2 (parsing/merging/sorting) time: " <<
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - start_time).count() <<
            " milliseconds\n";
    }
    

    std::string Bzhminmer_file = Bzhminmer_tmp_file + "." + std::to_string(nb_docs) + "c";
    std::string Bzhminmer_file_sparse = Bzhminmer_tmp_file + "_sparse." + std::to_string(nb_docs) + "c";

    // creating input reading buffer
    uint64_t elem_words = 1 + ((nb_docs + 63)/64); // number of uint64_t for 1 minmer + colors
    uint64_t buffer_bytes = 2 * constants::MB;
    uint64_t input_buffer_elems = buffer_bytes / (elem_words * sizeof(uint64_t)); // number of elems for a 2MB buffer 
    std::vector<uint64_t> input_buffer(input_buffer_elems * elem_words);

    ankerl::unordered_dense::set<uint64_t> unique_minmers;

    //DENSE FILE
    FILE* finput = fopen(Bzhminmer_file.c_str(), "rb");
    if (!finput) throw std::runtime_error("Cannot open parsing-result file");
    while (true) {
        size_t got = fread(input_buffer.data(), elem_words * sizeof(uint64_t), input_buffer_elems, finput);
        if (got == 0) break;
        for (size_t i = 0; i < got; ++i) {
            uint64_t minimizer = input_buffer[i * elem_words];
            unique_minmers.insert(minimizer);
        }
    }
    rewind(finput);

    uint64_t nb_elems_dense = unique_minmers.size();

    //SPARSE FILE (if exists)
    uint64_t shift = utils::sparse_colors_bits(nb_docs);
    uint64_t list_size_shift = 64 - shift;
    uint64_t granularity = 64 / utils::sparse_colors_bits(nb_docs);
    uint64_t granularity_log2 = utils::bits_needed(granularity) - 1; //avoid divisions later
    uint64_t buf_idx = 0;
    FILE* finput_sparse = fopen(Bzhminmer_file_sparse.c_str(), "rb");
    if (finput_sparse) {
        size_t got = fread(input_buffer.data(), sizeof(uint64_t), input_buffer_elems*elem_words, finput_sparse);
        while (buf_idx < got) {
            uint64_t minimizer = input_buffer[buf_idx++];
            unique_minmers.insert(minimizer);
            if (buf_idx >= got) {
                buf_idx = 0;
                got = fread(input_buffer.data(), sizeof(uint64_t), input_buffer_elems*elem_words, finput_sparse);
            }
            uint64_t list_size = input_buffer[buf_idx] >> list_size_shift;
            buf_idx += (list_size + granularity) / granularity;
            if (buf_idx >= got) {
                buf_idx = 0;
                got = fread(input_buffer.data(), sizeof(uint64_t), input_buffer_elems*elem_words, finput_sparse);
            }
        }

        rewind(finput_sparse);
    } 
    
    uint64_t nb_elems = unique_minmers.size();

    //STEP 3 : BUILDING MPHF ===================================================
    {
        if (build_parameters.verbose >= 1) std::cout << "[I] Step 3: building the MPHF for " << unique_minmers.size() << " minimizers\n";
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

        if (build_parameters.verbose >= 1) {
            std::cout << "[I] Step 3 (MPHF) time: " << 
                std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::high_resolution_clock::now() - start_time).count() << 
                " milliseconds\n";
        }
        
        if (build_parameters.verbose >= 2) {
            std::cout << "[II] MPHF Size : " << hf.num_bits()/8 << " Bytes = " << hf.num_bits()/constants::MB << " MB\n"; 
        }

        ankerl::unordered_dense::set<uint64_t>().swap(unique_minmers); //free memory
    }

    // common variables for CM & CS
    uint64_t ram_budget = static_cast<uint64_t>(build_parameters.max_ram_MB * 0.75); //75% of allocated RAM for CM & CS

    std::string minmer_to_cid_tmp_file = build_parameters.output_dirname + "/tmp/minmer_cid.bin";
    uint64_t output_buffer_elems = buffer_bytes / (2 * sizeof(uint64_t)); // number of (minmer-cid) for a 2MB buffer 

    uint64_t cid = -1; // starts at -1 so first unique set becomes 0


    //STEP 4.1 : COLORSETS DEDUP + STORING ====================================
    {
        if (build_parameters.verbose >= 1) std::cout << "[I] Step 4.1: colors deduplication and storing: ColorSets\n";
        start_time = std::chrono::high_resolution_clock::now();
        
        FILE* fout = fopen(minmer_to_cid_tmp_file.c_str(), "wb+");
        if (!fout) throw std::runtime_error("Cannot open tmp minmer-cid file");

        std::vector<uint64_t> output_buffer(output_buffer_elems * 2);

        uint64_t nb_words_colors = elem_words - 1;
        uint64_t* last_color_words = new uint64_t[nb_words_colors]();

        colorsets::builder cs_builder(
            build_parameters.output_dirname + "/colorsets", 
            nb_docs, 
            ram_budget * constants::MB,
            build_parameters.verbose
        );
        uint64_t output_pos = 0; // how many uint64_t currently in output buffer
        uint64_t processed = 0;  // number of elements processed so far
        
        // DENSE FILE
        while (processed < nb_elems_dense){
            // Read up to input_buf_elems elements
            size_t to_read = std::min(input_buffer_elems, nb_elems_dense - processed);
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

        buf_idx = 0;
        if (finput_sparse) {
            size_t got = fread(input_buffer.data(), sizeof(uint64_t), input_buffer_elems*elem_words, finput_sparse);

            // Use a vector to hold the last color set, resized as needed
            std::vector<uint64_t> last_sparse_colors;

            while (buf_idx < got) {
                uint64_t key = input_buffer[buf_idx++];
                if (buf_idx >= got) {
                    // current elem is split between minmer and colors
                    buf_idx = 0;
                    got = fread(input_buffer.data(), sizeof(uint64_t), input_buffer_elems*elem_words, finput_sparse);
                    if (got == 0) {
                        throw std::runtime_error("unexpected EOF (post minmer)");
                    }
                }

                uint64_t list_size = input_buffer[buf_idx] >> list_size_shift;
                nb_words_colors = (list_size + granularity) / granularity;

                if (buf_idx + nb_words_colors >= got) {
                    // Current color is split between buffers, preserve tail and refill
                    uint64_t tail_words = got - buf_idx;
                    memcpy(input_buffer.data(), &input_buffer[buf_idx], tail_words * sizeof(uint64_t));
                    buf_idx = 0;
                    size_t read_more = fread(input_buffer.data() + tail_words,
                                            sizeof(uint64_t),
                                            input_buffer_elems*elem_words - tail_words,
                                            finput_sparse);
                    got = tail_words + read_more;
                    if (got == 0) {
                        throw std::runtime_error("unexpected EOF (in colors)");
                    }
                }

                const uint64_t* color_words = &input_buffer[buf_idx];

                // Resize last_sparse_colors if needed
                if (last_sparse_colors.size() != nb_words_colors) {
                    last_sparse_colors.resize(nb_words_colors);
                }

                // Compare with last color set safely
                if (std::memcmp(color_words, last_sparse_colors.data(), nb_words_colors * sizeof(uint64_t)) != 0) {
                    std::vector<uint32_t> colors;
                    for (size_t c = 1; c <= list_size; c++) {
                        size_t word_idx = c >> granularity_log2; //= c / granularity;
                        size_t pos_in_word = c & (granularity - 1); //= c % granularity;
                        uint64_t word = color_words[word_idx];
                        
                        uint32_t color = static_cast<uint32_t>(
                            (word >> ((granularity - 1 - pos_in_word) * shift)) & ((1ULL << shift) - 1)
                        );
                        colors.push_back(color);
                    }
                    cs_builder.add_color_set(colors.data(), list_size);
                    cid++;
                    std::memcpy(last_sparse_colors.data(), color_words, nb_words_colors * sizeof(uint64_t));
                }

                // Output element: [key, cid+check_bits]
                uint64_t minmer_last_bits = key & ((1UL << build_parameters.b)-1);
                output_buffer[output_pos++] = hf(key); //hash of minmer
                output_buffer[output_pos++] = (cid << build_parameters.b) | minmer_last_bits;

                if (output_pos == output_buffer.size()) {
                    fwrite(output_buffer.data(), sizeof(uint64_t), output_pos, fout);
                    output_pos = 0;
                }

                buf_idx += nb_words_colors;
                processed++;
            }
        }

        // Flush remaining output
        if (output_pos > 0) {
            fwrite(output_buffer.data(), sizeof(uint64_t), output_pos, fout);
        }

        cs_builder.build(build_parameters.verbose);

        delete[] last_color_words;
        fclose(finput);
        fflush(fout);
        fclose(fout);
        
        if (build_parameters.verbose >= 1){
            std::cout << "[I] Step 4.1 (ColorSets) time: " << 
                std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::high_resolution_clock::now() - start_time).count() << 
                " milliseconds\n";
        }
    }

    //STEP 5 : COLORMAPPER  =============================================
    {
        if (build_parameters.verbose >= 1) std::cout << "[I] Step 4.2: mapping hash -> minimizer: ColorMapper\n";
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

        cm_builder.build(build_parameters.verbose);
        
        fclose(ffinal);

        if (build_parameters.verbose >= 1){
            std::cout << "[I] Step 4.2 (ColorMapper) time: " << 
                std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::high_resolution_clock::now() - start_time).count() << 
                " milliseconds\n";
        }
        
    }

    std::remove(Bzhminmer_file.c_str());
    std::remove(minmer_to_cid_tmp_file.c_str());
}

} // namespace minimizer
} // namespace kaminari


#endif // INDEX_BUILD_HPP