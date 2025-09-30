#ifndef KAMINARI_CODES_HPP
#define KAMINARI_CODES_HPP

#include "ext_bit_vector.hpp"
#include "ext_bit_parser.hpp"
#include "../bundled/biolib/include/bit_operations.hpp"

static inline void unary_encoder(ext_bit_vector& out, std::size_t x)
{
    while(x >= 64) { //unlikely
        out.push_back(uint64_t(0), 64);
        x -= 64;
    }
    uint64_t last = uint64_t(1) << x;
    out.push_back(last, x + 1);
}

static inline void gamma_encoder(ext_bit_vector& out, std::size_t x)
{
    auto xx = x + 1;
    std::size_t b = bit::msbll(static_cast<uint64_t>(xx));
    if (b >= 64) throw std::overflow_error("[gamma_encoder] Unable to gamma encode");
    unary_encoder(out, b);
    auto mask = (static_cast<uint64_t>(1) << b) - 1;
    out.push_back(xx & mask, b);
}

static inline void delta_encoder(ext_bit_vector& out, std::size_t x)
{
    auto xx = x + 1;
    std::size_t b = bit::msbll(static_cast<uint64_t>(xx));
    if (b >= 64) throw std::overflow_error("[delta_encoder] Unable to gamma encode");
    gamma_encoder(out, b);
    auto mask = (static_cast<uint64_t>(1) << b) - 1;
    out.push_back(xx & mask, b);
}



static inline std::size_t unary_decoder(ext_bit_parser& parser)
{
    return parser.parse_0();
}

static inline std::size_t gamma_decoder(ext_bit_parser& parser)
{
    std::size_t b = unary_decoder(parser);
    return (static_cast<std::size_t>(parser.parse_fixed(b)) | (std::size_t(1) << b)) - 1;
}

static inline std::size_t delta_decoder(ext_bit_parser& parser)
{
    std::size_t b = gamma_decoder(parser);
    return (parser.parse_fixed(b) | (std::size_t(1) << b)) - 1;
}

#endif // KAMINARI_CODES_HPP