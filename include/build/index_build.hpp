#ifndef INDEX_BUILD_HPP
#define INDEX_BUILD_HPP

#include "../../bundled/Minimizers/lib/BreiZHMinimizer.hpp"
#include "io_helper.hpp"
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
        false, 
        build_parameters.keep_tmp_files, 
        build_parameters.keep_tmp_files 
    );  

    vector<std::string>().swap(build_parameters.input_filenames); //free memory

    if (build_parameters.verbose >= 1){
        std::cout << "[I] Step 1 & 2 (parsing/merging/sorting) time: " <<
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - start_time).count() <<
            " milliseconds\n";
    }
    

    std::string Bzhminmer_file = Bzhminmer_tmp_file + "." + std::to_string(nb_docs) + "c.lz4";
    std::string Bzhminmer_file_sparse = Bzhminmer_tmp_file + "_sparse." + std::to_string(nb_docs) + "c.lz4";

    // creating input reading buffer
    uint64_t elem_words = 1 + ((nb_docs + 63)/64); // number of uint64_t for 1 minmer + colors
    uint64_t buffer_bytes = 2 * constants::MB;
    uint64_t input_buffer_elems = buffer_bytes / (elem_words * sizeof(uint64_t)); // number of elems for a 2MB buffer 
    std::vector<uint64_t> input_buffer(input_buffer_elems * elem_words);

    std::string temp_flat_file = "temp_minimizers_flat.bin";
    uint64_t nb_elems = 0;
    uint64_t nb_elems_dense;

    //STEP 3.a : BUILDING MPHF, gather data ===================================================
    {
        std::cout << "[I] Parsing and flattening minimizers to " << temp_flat_file << "...\n";
        
        // Open output stream with a buffer
        std::ofstream out_stream(temp_flat_file, std::ios::binary);
        if (!out_stream) {
            std::cerr << "Error: Could not create temp file.\n";
            exit(1);
        }

        const size_t WRITE_BUF_SIZE = 1024 * 1024; // 1M items (8MB)
        std::vector<uint64_t> write_buf;
        write_buf.reserve(WRITE_BUF_SIZE);

        auto flush_buffer = [&]() {
            if (!write_buf.empty()) {
                out_stream.write(reinterpret_cast<const char*>(write_buf.data()), 
                                 write_buf.size() * sizeof(uint64_t));
                write_buf.clear();
            }
        };

        // 1. DENSE FILE
        stream_reader* dense_parser = stream_reader_library::allocate(Bzhminmer_file);
        size_t got_dense = 0;
        while (true) {
            size_t items_read = dense_parser->read(input_buffer.data(), 
                                                   elem_words * sizeof(uint64_t), 
                                                   got_dense); 
            if (items_read == 0) break;
            for (size_t i = 0; i < items_read; ++i) { 
                uint64_t minimizer = input_buffer[i * elem_words];
                
                write_buf.push_back(minimizer);
                if (write_buf.size() >= WRITE_BUF_SIZE) flush_buffer();
                nb_elems++;
            }
        }
        delete dense_parser;
        nb_elems_dense = nb_elems;

        // 2. SPARSE FILE
        kaminari::helpers::SparseFileParser sparse_parser(Bzhminmer_file_sparse);
        if (sparse_parser.is_open()) {
            uint64_t minimizer;
            uint64_t bits = utils::sparse_colors_bits(nb_docs);
            while (sparse_parser.next_minimizer_only(minimizer, bits)) {
                write_buf.push_back(minimizer);
                if (write_buf.size() >= WRITE_BUF_SIZE) flush_buffer();
                nb_elems++;
            }
        }

        flush_buffer(); // Write remaining items
        out_stream.close();
        // Cleanup temp file
        std::remove(temp_flat_file.c_str());
    }

    //STEP 3.b : BUILDING MPHF  =========================================================
    {
        if (build_parameters.verbose >= 1) std::cout << "[I] Step 3: building the MPHF for " << nb_elems << " minimizers\n";
        start_time = std::chrono::high_resolution_clock::now();
        
        // Create the iterator we defined above
        kaminari::helpers::BinaryFileIterator file_it(temp_flat_file);

        // Standard pthash boilerplate
        int backup = dup(1);
        int redirect = open("/dev/null", O_WRONLY);
        fflush(stdout);
        dup2(redirect, 1);
        close(redirect);

        // BUILD: Pass the iterator and the count we just calculated
        hf.build_in_external_memory(file_it, nb_elems, get_pthash_options(build_parameters));

        fflush(stdout);
        dup2(backup, 1);
        close(backup);

        if (build_parameters.verbose >= 1) {
            std::cout << "[I] Step 3 (MPHF) time: " << 
                std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::high_resolution_clock::now() - start_time).count() << 
                " milliseconds\n";
        }
        
        if (build_parameters.verbose >= 2) {
            std::cout << "[II] MPHF Size : " << hf.num_bits()/8 << " Bytes = " << hf.num_bits()/constants::MB << " MB\n"; 
        }
    }

    
    // common variables for CM & CS
    uint64_t ram_budget = static_cast<uint64_t>(build_parameters.max_ram_MB * 0.75); //75% of allocated RAM for CM & CS

    std::string minmer_to_cid_tmp_file = build_parameters.output_dirname + "/tmp/minmer_cid.bin.lz4";
    uint64_t output_buffer_elems = buffer_bytes / (2 * sizeof(uint64_t)); // number of (minmer-cid) for a 2MB buffer 

    uint64_t cid = -1; // starts at -1 so first unique set becomes 0


    //STEP 4.1 : COLORSETS DEDUP + STORING ====================================
    {
        if (build_parameters.verbose >= 1) std::cout << "[I] Step 4.1: colors deduplication and storing: ColorSets\n";
        start_time = std::chrono::high_resolution_clock::now();
        
        stream_writer* fout = stream_writer_library::allocate(minmer_to_cid_tmp_file);
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
        stream_reader* dense_parser_step4 = stream_reader_library::allocate(Bzhminmer_file);
        while (processed < nb_elems_dense){
            // Read up to input_buf_elems elements
            size_t to_read = std::min(input_buffer_elems, nb_elems_dense - processed);
            size_t got = dense_parser_step4->read(input_buffer.data(), elem_words * sizeof(uint64_t), to_read);
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
                    fout->write(output_buffer.data(), sizeof(uint64_t), output_pos);
                    output_pos = 0;
                }
                processed++;
            }
        }
        delete dense_parser_step4;

        // SPARSE FILE
        kaminari::helpers::SparseFileParser sparse_parser_step4(Bzhminmer_file_sparse);
        if (sparse_parser_step4.is_open()) {
            uint64_t minimizer;
            std::vector<uint64_t> payload;
            std::vector<uint64_t> last_sparse_colors; // Cache for deduplication

            // Bit-unpacking configuration
            uint64_t bits_per_color = utils::sparse_colors_bits(nb_docs);
            uint64_t shift = bits_per_color;
            uint64_t list_size_shift = 64 - shift;
            uint64_t granularity = 64 / shift;
            uint64_t granularity_log2 = utils::bits_needed(granularity) - 1;

            // Loop efficiently extracts (Minimizer + Full Payload)
            while (sparse_parser_step4.next_element(minimizer, payload, bits_per_color)) {
                
                // 1. Deduplication: Check if payload (header + colors) changed
                if (last_sparse_colors != payload) {
                    // Extract colors from the packed payload
                    // payload[0] is the header, payload[1...] are packed colors
                    uint64_t list_size = payload[0] >> list_size_shift;
                    
                    std::vector<uint32_t> colors;
                    colors.reserve(list_size);

                    for (size_t c = 1; c <= list_size; c++) {
                        size_t word_idx = c >> granularity_log2;     // c / granularity
                        size_t pos_in_word = c & (granularity - 1);  // c % granularity
                        uint64_t word = payload[word_idx];
                        
                        uint32_t color = static_cast<uint32_t>(
                            (word >> ((granularity - 1 - pos_in_word) * shift)) & ((1ULL << shift) - 1)
                        );
                        colors.push_back(color);
                    }

                    // Add unique set
                    cs_builder.add_color_set(colors.data(), list_size);
                    cid++;
                    
                    // Update cache (vector assignment handles deep copy)
                    last_sparse_colors = payload;
                }

                // 2. Output Writing: [Hash(minmer), CID + check_bits]
                uint64_t minmer_last_bits = minimizer & ((1UL << build_parameters.b) - 1);
                output_buffer[output_pos++] = hf(minimizer); 
                output_buffer[output_pos++] = (cid << build_parameters.b) | minmer_last_bits;

                if (output_pos == output_buffer.size()) {
                    fout->write(output_buffer.data(), sizeof(uint64_t), output_pos);
                    output_pos = 0;
                }

                processed++;
            }
        }

        // Flush remaining output from buffer (handles both dense and sparse leftovers)
        if (output_pos > 0) {
            fout->write(output_buffer.data(), sizeof(uint64_t), output_pos);
        }

        cs_builder.build(build_parameters.verbose);

        delete[] last_color_words;
        delete fout;
        
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
            stream_reader* ffinal = stream_reader_library::allocate(minmer_to_cid_tmp_file);
            size_t processed = 0;
            while (processed < nb_elems) {
                size_t to_read = std::min(output_buffer_elems, nb_elems - processed);
                size_t got = ffinal->read(final_buffer.data(), 2 * sizeof(uint64_t), to_read);

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

            delete ffinal;
            low = high;
            high += cm_builder.flush();
        }

        // Final pass for the last chunk
        stream_reader* ffinal = stream_reader_library::allocate(minmer_to_cid_tmp_file);
        size_t processed = 0;
        while (processed < nb_elems) {
            size_t to_read = std::min(output_buffer_elems, nb_elems - processed);
            size_t got = ffinal->read(final_buffer.data(), 2 * sizeof(uint64_t), to_read);
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
        delete ffinal;

        cm_builder.build(build_parameters.verbose);
        
        if (build_parameters.verbose >= 1){
            std::cout << "[I] Step 4.2 (ColorMapper) time: " << 
                std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::high_resolution_clock::now() - start_time).count() << 
                " milliseconds\n";
        }
        
    }

    //CLEANING TMP FILES =======================================================
    if (!build_parameters.keep_tmp_files) {
        std::remove(Bzhminmer_file.c_str());
        std::remove(Bzhminmer_file_sparse.c_str());
        std::remove(minmer_to_cid_tmp_file.c_str());
    }
}

} // namespace minimizer
} // namespace kaminari


#endif // INDEX_BUILD_HPP