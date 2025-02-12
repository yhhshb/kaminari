#ifndef KAMINARI_INDEX_HPP
#define KAMINARI_INDEX_HPP

#include <zlib.h>
extern "C" {
#include "../bundled/kseq.h"
}

#include <iostream>
#include <mutex>
#include <shared_mutex>
#include <condition_variable>
#include <deque>
#include <stack>
#include <sys/stat.h>


#include "../bundled/pthash/include/pthash.hpp"
#include "../bundled/biolib/bundled/prettyprint.hpp"
#include "../bundled/biolib/include/bit_vector.hpp"
#include "../bundled/biolib/include/elias_fano.hpp"
#include "../bundled/biolib/include/iterator/standalone_iterator.hpp"
#include "../bundled/biolib/include/iterator/sorted_merge_iterator.hpp"
#include "../bundled/unordered_dense/include/ankerl/unordered_dense.h"

#include "constants.hpp"
#include "minimizer.hpp"
#include "utils.hpp"
#include "build_options.hpp"
#include "query_options.hpp"
#include "compact_vector.hpp"

KSEQ_INIT(gzFile, gzread)

namespace kaminari {
namespace minimizer {

#define CLASS_HEADER template <class ColorClasses, class ColorMapper>
#define METHOD_HEADER index<ColorClasses, ColorMapper>

template <typename T>
struct scored {
    T item;
    uint32_t score;
};
typedef scored<uint32_t> scored_id;

struct element {
    uint64_t minmer;    // le minimizer
    uint64_t* colors; // les couleurs associées
    int n_blocks;
};







CLASS_HEADER
class index
{
    public:
        typedef uint64_t minimizer_t;
        typedef typename ColorClasses::color_t color_t;
        typedef typename kaminari::query::options_t options_t;

        index();
        index(const build::options_t& build_parameters);

        //following methods are explicitly instantiated in src/psa/
        std::vector<color_t> query_union_threshold(char const * const q, const std::size_t l, options_t& opts) const noexcept;
        std::vector<scored_id> ranking_query_union_threshold(char const * const q, const std::size_t l, options_t& opts) const noexcept;
        
        void memory_breakdown(std::ostream& out) const noexcept;
        
        template <class Visitor>
        void visit(Visitor& visitor);

    private:
        typedef pthash::build_configuration pthash_opt_t;
        //typedef pthash::phobic<pthash::xxhash128> pthash_minimizers_mphf_t; //TODO: not sure about visit yet
        typedef pthash::dense_partitioned_phf<  
            pthash::xxhash128, //murmurhash2_64                     
            pthash::opt_bucketer,               
            pthash::mono_EF,                    
            true,
            pthash::pthash_search_type::add_displacement  
            >
            pthash_minimizers_mphf_t;

        template <class Iterator>
        class pthash_input_iterator {
            public:
                pthash_input_iterator(Iterator mm_itr) : m_iterator(mm_itr) {}
                void operator++() {++m_iterator;}
                typename Iterator::value_type operator*() const {return (*m_iterator);}
            private:
                Iterator m_iterator;
        };
        pthash_opt_t get_pthash_options(const build::options_t& build_parameters);

        
        uint64_t get_file_size(const std::string& filen) const;
        static bool Sgreater_func (const element& e1, const element& e2);
        static bool Sequal_func (const element& e1, const element& e2);
        void process(const std::string& ifile, std::vector<element>& buffer, int n_colors);

        void build(const build::options_t& build_parameters);



        //following methods are explicitly instantiated in src/psa/files
        //with colorsclasses being from hybrid.hpp and color mapper being pthash::compact_vector
        std::vector<scored_id> ranking_dense_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, color_t>>&& color_id_itrs, uint64_t threshold) const noexcept;
        std::vector<scored_id> ranking_mixed_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, color_t>>&& color_id_itrs, uint64_t threshold) const noexcept;

        std::vector<color_t> union_dense_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, color_t>>&& color_id_itrs, uint64_t threshold) const noexcept;
        std::vector<color_t> union_mixed_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, color_t>>&& color_id_itrs, uint64_t threshold) const noexcept;
        
        std::vector<std::string> m_filenames;
        uint8_t k;
        uint8_t m;
        uint8_t b;
        uint64_t seed;
        bool canonical;
        double pthash_constant;
        ColorClasses m_ccs; // colors
        ColorMapper m_map; // map between mphf values and color classes
        pthash_minimizers_mphf_t hf; // minimizer mphf
};

CLASS_HEADER
METHOD_HEADER::index()
    :
    k(0),
    m(0),
    b(0),
    seed(0),
    canonical(false),
    pthash_constant(0)
{}

CLASS_HEADER
METHOD_HEADER::index(const build::options_t& build_parameters)
    : 
    m_filenames(build_parameters.input_filenames),
    k(build_parameters.k),
    m(build_parameters.m),
    b(build_parameters.b),
    seed(build_parameters.seed),
    canonical(build_parameters.canonical),
    pthash_constant(build_parameters.pthash_constant)
{
    build(build_parameters);
}

CLASS_HEADER
void 
METHOD_HEADER::memory_breakdown(std::ostream& out) const noexcept
{
    libra scale;
    scale.visit(m_filenames);
    out << "The list of input filenames weights: " << scale.get_byte_size() << " Bytes\n";
    scale.reset();
    out << "The MPHF of minimizers weights: " << hf.num_bits() / 8 << " Bytes\n";
    scale.visit(m_ccs);
    out << "Colors weight: " << scale.get_byte_size() << " Bytes\n";
    scale.reset();
    //TODO:fix this
    //scale.visit(m_map);
    out << "The mapping from minimizers to colors weights: " << m_map.num_bytes() << " Bytes\n";
    //scale.reset();
}

CLASS_HEADER
template <class Visitor>
void
METHOD_HEADER::visit(Visitor& visitor)
{
    visitor.visit(m_filenames);
    visitor.visit(k);
    visitor.visit(m);
    visitor.visit(b);
    visitor.visit(seed);
    visitor.visit(canonical);
    visitor.visit(pthash_constant);
    visitor.visit(hf); // lphash mphf
    visitor.visit(m_ccs); // colors
    visitor.visit(m_map); // map between mphf values and color classes
}

CLASS_HEADER
typename METHOD_HEADER::pthash_opt_t
METHOD_HEADER::get_pthash_options(const build::options_t& build_parameters)
{
    pthash_opt_t opts;
    opts.seed = build_parameters.seed;
    opts.lambda = build_parameters.pthash_constant; // (too slow = try decreasing), higher lambda : more space efficient 
    opts.alpha = 0.97; //was 0.94
    opts.search = pthash::pthash_search_type::add_displacement;
    opts.avg_partition_size = 3000;
    opts.verbose = (build_parameters.verbose > 0);
    
    opts.ram = build_parameters.max_ram * constants::GB;
    opts.num_threads = build_parameters.nthreads;
    opts.tmp_dir = build_parameters.tmp_dir;

    opts.dense_partitioning = true;

    return opts;
}

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
void METHOD_HEADER::process(const std::string& ifile, std::vector<element>& buffer, int n_colors)
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

    printf("(II) file size in bytes  : %llu\n", size_bytes);
    printf("(II) # uint64_t elements : %llu\n", n_elements);
    printf("(II) # minimizers        : %llu\n", n_minimizr);
    printf("(II) # of colors         : %d\n",   n_colors  );
    printf("(II) # uint64_t/colors   : %llu\n", n_uint64_c);

    ////////////////////////////////////////////////////////////////////////////////////

    FILE* fi = fopen( ifile.c_str(), "r" );
    if( fi == NULL )
    {
        printf("(EE) An error corrured while openning the file (%s)\n", ifile.c_str());
        printf("(EE) Error location : %s %d\n", __FILE__, __LINE__);
        exit( EXIT_FAILURE );
    }

    buffer.resize(n_minimizr);
    for (int i = 0; i < n_minimizr; i += 1){
        buffer[i].colors = new uint64_t[n_uint64_c];
        const int n_minmer_reads = fread(&buffer[i].minmer, sizeof(uint64_t), 1, fi);
        const int n_colors_reads = fread(buffer[i].colors, sizeof(uint64_t), n_uint64_c, fi);
        buffer[i].n_blocks = n_uint64_c;
    }

    
    fclose( fi );

    std::sort( buffer.begin(), buffer.end(), &Sgreater_func);
}

CLASS_HEADER
void
METHOD_HEADER::build(const build::options_t& build_parameters)
{
    //STEP 1 : PARSE FILES =====================================================
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::string sorted_file = "/home/vlevallo/tmp/test_bertrand/result.4c";

    std::vector<element> sorted_buffer;
    process(sorted_file, sorted_buffer, m_filenames.size());

    std::vector<minimizer_t> unique_minmers(sorted_buffer.size());
    for (int i = 0; i < sorted_buffer.size(); i += 1){
        unique_minmers[i] = sorted_buffer[i].minmer;
    }

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
        // TODO: ceil(log2(hf.num_keys())) depends on the number of unique minmer, should depend on the number of distinct colors instead, but should not bug because nb_distinct_colors <= nb_unique_minmers

        color_t cid = 0;
        color_t cid_with_parity = 0;

        uint64_t minimizer;
        uint64_t mp_idx;


        for (int i = 0; i < sorted_buffer.size(); i += 1){
            //add first time seen color to ColorClasses
            //was under binary form, need to translate it
            auto binary_color = sorted_buffer[i].colors;
            std::vector<color_t> color_list;  
            for (std::size_t block_idx = 0; block_idx < sorted_buffer[i].n_blocks; ++block_idx) {
                uint64_t block = binary_color[block_idx];
                for (std::size_t bit_pos = 0; bit_pos < 64; ++bit_pos) {
                    if (block & (1ULL << bit_pos)) { 
                        color_list.push_back(block_idx * 64 + bit_pos);
                    }
                }
            }
            cbuild.add_color_set(color_list.data(), color_list.size());

            //link minmer to cid
            minimizer = sorted_buffer[i].minmer;
            mp_idx = hf(minimizer);
            cid_with_parity = (cid << build_parameters.b) | ( minimizer & ((1UL << build_parameters.b)-1) );
            m_map_builder.set(mp_idx, cid_with_parity);

            while (i+1 < sorted_buffer.size() and Sequal_func(sorted_buffer[i], sorted_buffer[i+1])){
                //ignore redundant color
                i += 1;
                //link minmer to cid
                minimizer = sorted_buffer[i].minmer;
                mp_idx = hf(minimizer);
                cid_with_parity = (cid << build_parameters.b) | ( minimizer & ((1UL << build_parameters.b)-1) );
                m_map_builder.set(mp_idx, cid_with_parity);
            }

            cid += 1;
        }

        for (int i = 0; i < sorted_buffer.size(); i++) {
            delete[] sorted_buffer[i].colors; 
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

#endif // KAMINARI_INDEX_HPP