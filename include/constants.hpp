#ifndef KAMINARI_CONSTANTS_HPP
#define KAMINARI_CONSTANTS_HPP

#include <cstddef>
#include <zlib.h>
extern "C" {
#include "../bundled/kseq.h"
}
#include "compile_constants.tdp"

#include "../bundled/biolib/include/bit_vector.hpp"
#include "../bundled/biolib/include/packed_vector.hpp"
#include "../bundled/biolib/include/bit_parser.hpp"
#include "../bundled/biolib/include/elias_fano.hpp"
#include "../bundled/biolib/include/codes.hpp"
#include "../bundled/biolib/include/external_memory_vector.hpp"
#include "../bundled/biolib/include/counting_iterator.hpp"
#include "../bundled/biolib/include/ordered_unique_sampler.hpp"
#include "../bundled/biolib/include/member_iterator.hpp"
#include "../bundled/biolib/include/io.hpp"
#include "../bundled/biolib/include/logtools.hpp"
#include "../bundled/biolib/include/hash.hpp"
#include "../bundled/pthash/include/pthash.hpp"
// #include "../bundled/cppitertools/groupby.hpp"

KSEQ_INIT(gzFile, gzread)

namespace kaminari {

typedef bit::vector<uint64_t> bit_vector;
typedef bit::packed::vector<uint64_t> packed_vector;
typedef bit::rs::array<bit_vector, 64, 8, false, false> ranked_bit_vector;
typedef bit::ef::array ef_sequence;
typedef bit::parser<uint64_t> bit_parser;
typedef emem::external_memory_vector<uint32_t, false> emem_colors;
typedef iterators::counting_iterator<uint32_t> dummy_itr_t;
typedef io::mut_saver saver;
typedef io::loader loader;
typedef logging_tools::libra libra;
typedef hash::hash64 hash64;

typedef pthash::build_configuration pthash_opt_t;
typedef pthash::single_phf<pthash::murmurhash2_64, pthash::dictionary_dictionary, true> pthash_minimizers_mphf_t;

typedef uint64_t minimizer_t;
typedef uint32_t color_t;

struct opt_t {
    using fn_t = std::vector<std::string>;
    fn_t input_filenames;
    std::string output_filename;
    std::string tmp_dir;
    uint8_t k;
    uint8_t m;
    uint8_t nthreads;
    std::size_t max_ram;
    uint64_t seed;
    double pthash_constant;
    bool canonical;
    bool check;
    std::size_t verbose;
};

namespace constants {

static const std::size_t MAX_KMER_SIZE = sizeof(kmer_t) * 4;
static const std::size_t GB = 1000 * 1000 * 1000;

}

namespace util {

std::string get_tmp_filename(const std::string& tmp_dirname, const std::string& prefix, uint64_t run_identifier);

} // namespace util

namespace non_standard {

template <class Iterator>
class pthash_input_iterator {
public:
    pthash_input_iterator(Iterator mm_itr) : m_iterator(mm_itr) {}

    void operator++() {++m_iterator;}

    typename Iterator::value_type operator*() const {return (*m_iterator);}

private:
    Iterator m_iterator;
};

} // namespace non_standard

} // namespace kaminari

#endif // KAMINARI_CONSTANTS_HPP
