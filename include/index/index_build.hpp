#ifndef INDEX_BUILD_HPP
#define INDEX_BUILD_HPP

#include <cstdlib>  // For system()

namespace kaminari {
namespace minimizer {

#define CLASS_HEADER template <class ColorClasses, class ColorMapper>
#define METHOD_HEADER index<ColorClasses, ColorMapper>

CLASS_HEADER
uint64_t METHOD_HEADER::get_file_size(const std::string& filen) const {
    struct stat file_status;
    if (stat(filen.c_str(), &file_status) < 0) {
        return -1;
    }
    return file_status.st_size;
}

CLASS_HEADER
bool METHOD_HEADER::Sgreater_func (const element& e1, const element& e2)
{
    for(int i = 0; i < e1.n_blocks; i += 1){
        if( e1.colors[i] != e2.colors[i] ){
            return e1.colors[i] < e2.colors[i];
        }
    }
    return e1.minmer > e2.minmer;
}

CLASS_HEADER
bool METHOD_HEADER::Sequal_func (const element& e1, const element& e2)
{
    for(int i = 0; i < e1.n_blocks; i += 1){
        if( e1.colors[i] != e2.colors[i] ) return false;
    }
    return true;
}

CLASS_HEADER
void METHOD_HEADER::process(const std::string& ifile, std::vector<element>& min_col, int n_colors)
{
    const uint64_t size_bytes  = get_file_size(ifile);
    const uint64_t n_elements  = size_bytes / sizeof(uint64_t);
    uint64_t n_minimizr;
    uint64_t n_uint64_c;

    if( n_colors < 64  ){
        n_uint64_c  = 1;
    }else{
        n_uint64_c  = ((n_colors + 63) / 64);
    }

    n_minimizr  = n_elements / (1 + n_uint64_c);
    min_col.resize(n_minimizr);

    printf("(II) file size in bytes  : %llu\n", size_bytes);
    printf("(II) # uint64_t elements : %llu\n", n_elements);
    printf("(II) # minimizers        : %llu\n", n_minimizr);
    printf("(II) # of colors         : %d\n",   n_colors  );
    printf("(II) # uint64_t/colors   : %llu\n", n_uint64_c);

    ////////////////////////////////////////////////////////////////////////////////////

    uint64_t buffer[4096];

    FILE* fi = fopen( ifile.c_str(), "r" );
    if( fi == NULL )
    {
        printf("(EE) An error corrured while openning the file (%s)\n", ifile.c_str());
        printf("(EE) Error location : %s %d\n", __FILE__, __LINE__);
        exit( EXIT_FAILURE );
    }

    int cnt = 0;
    uint64_t n_data = fread(buffer, sizeof(uint64_t), 4096, fi);
    
    for (int i = 0; i < n_minimizr; i += 1){
        min_col[i].minmer = buffer[cnt];
        cnt += 1;
        if (cnt == 4096){
            n_data = fread(buffer, sizeof(uint64_t), 4096, fi);
            cnt = 0;
        }

        min_col[i].colors = new uint64_t[n_uint64_c];
        for (int j = 0; j < n_uint64_c; j += 1){
            min_col[i].colors[j] = buffer[cnt];
            cnt += 1;
            if (cnt == 4096){
                n_data = fread(buffer, sizeof(uint64_t), 4096, fi);
                cnt = 0;
            }
        }

        min_col[i].n_blocks = n_uint64_c;
    }

    fclose( fi );
    
    std::sort( min_col.begin(), min_col.end(), &Sgreater_func);
}

CLASS_HEADER
void
METHOD_HEADER::build(const build::options_t& build_parameters)
{
    //STEP 1 : PARSE FILES =====================================================
    if (build_parameters.verbose > 0) std::cerr << "Step 1: parsing " << build_parameters.input_filenames.size() << " files and merging results\n";

    auto start_time = std::chrono::high_resolution_clock::now();

    auto allocated = get_file_size(m_filenames[0]);
    allocated = (allocated/1024/1024 < 10) ? 16 : 1024; // 16MB for bacteria/small genomes, 1GB for other genomes

    if (build_parameters.fof_filename != "") {
        std::string build_dir = utils::getExecutablePath();
        std::string command = 
            build_dir + 
            "/BreiZHMinimizer" + 
            " -f " + build_parameters.fof_filename +
            " -t " + std::to_string(build_parameters.nthreads) +
            " --MB " + std::to_string(allocated);
        system(command.c_str());
    } else {
        std::cerr << "(EE) Please use a file of files if --metagenome is not used\n";
        exit(EXIT_FAILURE);
    }
    
    std::vector<element> sorted_min_cols;
    std::string Bzhminmer_file = "result." + std::to_string(build_parameters.input_filenames.size()) + "c";
    process(Bzhminmer_file, sorted_min_cols, m_filenames.size());

    ankerl::unordered_dense::set<minimizer_t> unique_minmers;
    for (int i = 0; i < sorted_min_cols.size(); i += 1){
        unique_minmers.insert(sorted_min_cols[i].minmer);
    }

    std::cerr << "Number of unique minimizers after sorting and before MPHF: " << unique_minmers.size() << std::endl;

    std::cerr << "Time for reading colors + sort them: " 
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - start_time
                    ).count() 
                << " milliseconds\n";

    //STEP 3 : BUILDING MPHF ===================================================
    {
        if (build_parameters.verbose > 0) std::cerr << "Step 3: building the MPHF for " << unique_minmers.size() << " minimizers\n";
        start_time = std::chrono::high_resolution_clock::now();
        //auto pt_itr = pthash_input_iterator<decltype(unique_minimizers)::const_iterator>(unique_minimizers.cbegin());
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

        std::cerr << "Time for MPHF build: " 
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - start_time
                    ).count() 
                << " milliseconds\n";

        std::cerr << "size of MPHF : " << hf.num_bits()/8 << " Bytes\n";
    }

    //STEP 4 : LIST DEDUPLICATION + MAPPING ====================================
    {
        if (build_parameters.verbose > 0) std::cerr << "Step 4: list deduplication and mapping\n";
        start_time = std::chrono::high_resolution_clock::now();
        
        typename ColorClasses::builder cbuild(m_filenames.size(), build_parameters.verbose);
        kaminari::compact_vector::builder m_map_builder(hf.num_keys(), ceil(log2(hf.num_keys()))+build_parameters.b); //1bit for check

        color_t cid = 0;
        color_t cid_with_parity = 0;

        uint64_t minimizer;
        uint64_t mp_idx;


        for (int i = 0; i < sorted_min_cols.size(); i += 1){
            //add first time seen color to ColorClasses
            //was under binary form, need to translate it
            auto binary_color = sorted_min_cols[i].colors;
            std::vector<color_t> color_list;  
            for (std::size_t block_idx = 0; block_idx < sorted_min_cols[i].n_blocks; ++block_idx) {
                uint64_t block = binary_color[block_idx];
                for (std::size_t bit_pos = 0; bit_pos < 64; ++bit_pos) {
                    if (block & (1ULL << bit_pos)) { 
                        color_list.push_back(block_idx * 64 + bit_pos);
                    }
                }
            }
            cbuild.add_color_set(color_list.data(), color_list.size());

            //link minmer to cid
            minimizer = sorted_min_cols[i].minmer;
            mp_idx = hf(minimizer);
            cid_with_parity = (cid << build_parameters.b) | ( minimizer & ((1UL << build_parameters.b)-1) );
            m_map_builder.set(mp_idx, cid_with_parity);

            while (i+1 < sorted_min_cols.size() and Sequal_func(sorted_min_cols[i], sorted_min_cols[i+1])){
                //ignore redundant color
                i += 1;
                //link minmer to cid
                minimizer = sorted_min_cols[i].minmer;
                mp_idx = hf(minimizer);
                cid_with_parity = (cid << build_parameters.b) | ( minimizer & ((1UL << build_parameters.b)-1) );
                m_map_builder.set(mp_idx, cid_with_parity);
            }

            cid += 1;
        }

        for (int i = 0; i < sorted_min_cols.size(); i++) {
            delete[] sorted_min_cols[i].colors; 
        }

        //shrink color mapper because we store cids so we only need log2(cid) bits
        m_map_builder.shrink(ceil(log2(hf.num_keys())) - ceil(log2(cid)));

        cbuild.build(m_ccs);
        m_map_builder.build(m_map);
        //compact_vector m_map where m_map[ hf(minmer) ] = color_id 
    }

    std::cerr << "Time for dedup + mapping: " 
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - start_time
                    ).count() 
                << " milliseconds\n";

    if (build_parameters.verbose > 0) {
        std::cerr << "Number of ids: " << m_ccs.num_docs() << "\n";
        std::cerr << "Number of colors (lists of ids):" << m_ccs.num_color_classes() << "\n";
    } 
}


#undef CLASS_HEADER
#undef METHOD_HEADER


} // namespace minimizer
} // namespace kaminari


#endif // INDEX_BUILD_HPP