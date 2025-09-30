#ifndef KAMINARI_INDEX_HPP
#define KAMINARI_INDEX_HPP

#include <zlib.h>
extern "C" {
#include "../../bundled/kseq.h"
}

#include <iostream>
#include <mutex>
#include <shared_mutex>
#include <condition_variable>
#include <deque>
#include <stack>
#include <sys/stat.h>

#include "../../bundled/pthash/include/pthash.hpp"
#include "../../bundled/biolib/bundled/prettyprint.hpp"
#include "../../bundled/biolib/include/bit_vector.hpp"
#include "../../bundled/biolib/include/elias_fano.hpp"
#include "../../bundled/biolib/include/iterator/standalone_iterator.hpp"
#include "../../bundled/biolib/include/iterator/sorted_merge_iterator.hpp"
#include "../../bundled/unordered_dense/include/ankerl/unordered_dense.h"

#include "../constants.hpp"
#include "../minimizer.hpp"
#include "../utils.hpp"
#include "../build_options.hpp"
#include "../query_options.hpp"
#include "../compact_vector.hpp"

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

struct element { //for build()
    uint64_t minmer;    // le minimizer
    uint64_t* colors; // les couleurs associées
    size_t n_blocks;
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

        
        ////
        // in index_build.hpp, classic index with Minimizer library
        ////
        uint64_t get_file_size(const std::string& filen) const;
        static bool Sgreater_func (const element& e1, const element& e2);
        static bool Sequal_func (const element& e1, const element& e2);
        void process(const std::string& ifile, std::vector<element>& buffer, int n_colors);

        void build(const build::options_t& build_parameters);


        ////
        // in psa/ directory, query the index
        //following methods are explicitly instantiated in include/psa/files
        //with colorsclasses being from hybrid.hpp and color mapper being pthash::compact_vector
        std::vector<scored_id> ranking_dense_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, color_t>>&& color_id_itrs, uint64_t threshold) const noexcept;
        std::vector<scored_id> ranking_mixed_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, color_t>>&& color_id_itrs, uint64_t threshold) const noexcept;
        
        ////
        // members
        ////
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



#undef CLASS_HEADER
#undef METHOD_HEADER

} // namespace minimizer
} // namespace kaminari

#endif // KAMINARI_INDEX_HPP