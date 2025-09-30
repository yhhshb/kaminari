#ifndef KAMINARI_INDEX_HPP
#define KAMINARI_INDEX_HPP

#include <zlib.h>
extern "C" {
#include "../../bundled/kseq.h"
}

#include <iostream>

#include "../bundled/pthash/include/pthash.hpp"
#include "../bundled/biolib/bundled/prettyprint.hpp"
#include "../bundled/biolib/include/bit_vector.hpp"
#include "../bundled/biolib/include/elias_fano.hpp"
#include "../bundled/biolib/include/iterator/standalone_iterator.hpp"
#include "../bundled/biolib/include/iterator/sorted_merge_iterator.hpp"
#include "../bundled/unordered_dense/include/ankerl/unordered_dense.h"

#include "constants.hpp"
#include "query/minimizer.hpp"
#include "build/build_options.hpp"
#include "query/query_options.hpp"
#include "colormapper.hpp"
#include "query/colorsets_accessor.hpp"

KSEQ_INIT(gzFile, gzread)

namespace kaminari {
namespace minimizer {

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

class index
{
    public:
        typedef typename kaminari::query::options_t options_t;

        index();
        index(build::options_t& build_parameters);

        //following methods are explicitly instantiated in src/psa/
        std::vector<scored_id> ranking_query_union_threshold(char const * const q, const std::size_t l, options_t& opts) const noexcept;
        
        void memory_breakdown(std::ostream& out) const noexcept;
        
        template <class Visitor>
        void visit(Visitor& visitor)
        {
            visitor.visit(nb_docs);
            visitor.visit(k);
            visitor.visit(m);
            visitor.visit(b);
            visitor.visit(seed);
            visitor.visit(canonical);
            visitor.visit(pthash_constant);
            visitor.visit(hf); // lphash mphf
            //visitor.visit(m_ccs); // colors
            //visitor.visit(m_map); // map between mphf values and color classes
        }

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
        pthash_opt_t get_pthash_options(build::options_t& build_parameters);

        
        ////
        // in index_build.hpp, classic index with Minimizer library
        ////
        uint64_t get_file_size(const std::string& filen) const;
        static bool Sgreater_func (const element& e1, const element& e2);
        static bool Sequal_func (const element& e1, const element& e2);
        void process(const std::string& ifile, std::vector<element>& buffer, int n_colors);

        void build(build::options_t& build_parameters);

        //std::vector<scored_id> ranking_mixed_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, color_t>>&& color_id_itrs, uint64_t threshold) const noexcept;
        
        ////
        // members
        ////
        uint64_t nb_docs;
        uint8_t k;
        uint8_t m;
        uint8_t b;
        uint64_t seed;
        bool canonical;
        double pthash_constant;
        colorsets_accessor m_colorsets; // colors
        colormapper m_colormapper; // map between mphf values and color classes
        pthash_minimizers_mphf_t hf; // minimizer mphf
};


} // namespace minimizer
} // namespace kaminari

#endif // KAMINARI_INDEX_HPP