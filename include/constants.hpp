#pragma once

#include <cstddef>
#include "compile_constants.tdp"
#include "../bundled/biolib/include/bit_vector.hpp"
#include "../bundled/biolib/include/bit_parser.hpp"
#include "../bundled/biolib/include/elias_fano.hpp"
#include "../bundled/biolib/include/codes.hpp"

namespace kaminari {

typedef bit::vector<uint64_t> bit_vector;
typedef bit::ef::array ef_sequence;
typedef bit::parser<uint64_t> bit_parser;

namespace constants {

static const std::size_t MAX_KMER_SIZE = sizeof(kmer_t) * 4;

}
}
