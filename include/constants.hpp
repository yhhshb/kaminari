#ifndef KAMINARI_CONSTANTS_HPP
#define KAMINARI_CONSTANTS_HPP

#include <cstdint>
#include <cstddef>
#include <zlib.h>
extern "C" {
#include "../bundled/kseq.h"
}
#include "compile_constants.tdp"
#include "../bundled/biolib/include/external_memory_vector.hpp"
#include "../bundled/biolib/include/ordered_unique_sampler.hpp"
#include "../bundled/biolib/include/iterator/member_iterator.hpp"
#include "../bundled/biolib/include/io.hpp"
#include "../bundled/biolib/include/logtools.hpp"
#include "../bundled/biolib/include/hash.hpp"

KSEQ_INIT(gzFile, gzread)

namespace kaminari {

// typedef bit::rs::array<bit_vector, 64, 8, false, false> ranked_bit_vector;
typedef emem::external_memory_vector<uint32_t, false> emem_colors;
typedef io::mut_saver saver;
typedef io::loader loader;
typedef logging_tools::libra libra;
typedef hash::hash64 hash64;

namespace constants {

static const std::size_t MAX_KMER_SIZE = sizeof(kmer_t) * 4;
static const std::size_t GB = 1000 * 1000 * 1000;

}

namespace util {

std::string get_tmp_filename(const std::string& tmp_dirname, const std::string& prefix, uint64_t run_identifier);

} // namespace util

} // namespace kaminari

#endif // KAMINARI_CONSTANTS_HPP
