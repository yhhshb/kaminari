#ifndef KAMINARI_CONSTANTS_HPP
#define KAMINARI_CONSTANTS_HPP

#include <cstddef>
#include "compile_constants.tdp"

#include "../bundled/biolib/include/bit_vector.hpp"
#include "../bundled/biolib/include/bit_parser.hpp"
#include "../bundled/biolib/include/elias_fano.hpp"
#include "../bundled/biolib/include/codes.hpp"
#include "../bundled/biolib/include/external_memory_vector.hpp"
#include "../bundled/biolib/include/counting_iterator.hpp"
#include "../bundled/biolib/include/io.hpp"
#include "../bundled/biolib/include/logtools.hpp"
#include "../bundled/lphash/lib/include/partitioned_mphf.hpp"
#include "../bundled/lphash/main/include/constants.hpp"

namespace kaminari {

typedef bit::vector<uint64_t> bit_vector;
typedef bit::rs::array<bit_vector, 64, 8, false> ranked_bit_vector;
typedef bit::ef::array ef_sequence;
typedef bit::parser<uint64_t> bit_parser;
typedef emem::external_memory_vector<uint32_t, false> emem_colors;
typedef iterators::counting_iterator<uint32_t> dummy_itr_t;
typedef io::mut_saver saver;
typedef io::loader loader;
typedef logging_tools::libra libra;

typedef lphash::mphf::partitioned lphash_mphf_t;
typedef lphash::mphf::interface::configuration lphash_configuration_t;

struct opt_t {
    using fn_t = std::vector<std::string>;
    fn_t input_filenames;
    std::string output_filename;
    std::string tmp_dir;
    uint8_t k;
    uint8_t m;
    uint8_t nthreads;
    std::size_t max_ram;
    bool check;
    bool verbose;
};

namespace constants {

static const std::size_t MAX_KMER_SIZE = sizeof(kmer_t) * 4;
static const auto lphash_c = lphash::constants::c;

}
}

#endif // KAMINARI_CONSTANTS_HPP
