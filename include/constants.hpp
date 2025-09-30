#ifndef KAMINARI_CONSTANTS_HPP
#define KAMINARI_CONSTANTS_HPP

#include <cstdint>
#include <cstddef>

#include "../bundled/biolib/include/external_memory_vector.hpp"
#include "../bundled/biolib/include/ordered_unique_sampler.hpp"
#include "../bundled/biolib/include/iterator/member_iterator.hpp"
#include "../bundled/biolib/include/io.hpp"
#include "../bundled/biolib/include/logtools.hpp"
#include "../bundled/biolib/include/hash.hpp"

namespace kaminari {

// typedef bit::rs::array<bit_vector, 64, 8, false, false> ranked_bit_vector;
typedef emem::external_memory_vector<std::pair<uint64_t, uint32_t>> emem_t;
typedef emem::external_memory_vector<std::pair<std::vector<uint32_t>, uint64_t>> colors_to_minmer;
typedef io::mut_saver saver;
typedef io::loader loader;
typedef logging_tools::libra libra;
typedef hash::hash64 hash64;
typedef hash::double_hash64 double_hash64;



namespace constants {

static const std::size_t MAX_KMER_SIZE = 32;

static const std::size_t KB = 1024;
static const std::size_t MB = 1024 * KB;
static const std::size_t GB = 1024 * MB;

} // namespace constants

} // namespace kaminari

#endif // KAMINARI_CONSTANTS_HPP
