#include "../include/ext_bit_parser.hpp" 

ext_bit_parser::ext_bit_parser() 
    : _data(nullptr), block_size(0), idx(0), buffer(0), available(0) 
{}

ext_bit_parser::ext_bit_parser(uint64_t const * const data, std::size_t size, std::size_t starting_bit_position)
    : _data(data), block_size(size)
{
    reset(starting_bit_position);
}

ext_bit_parser::ext_bit_parser(mymm::immap<uint64_t>& mapped, std::size_t starting_bit_position)
            : _data(mapped.data()), block_size(64)
        {
            reset(starting_bit_position);
        }

/* return the next l bits from the current position and advance by l bits */
uint64_t 
ext_bit_parser::parse_fixed(std::size_t l)
{
    uint64_t val = 0;
    if (l >= ::bit::size(val)) throw std::length_error("[bit::ext_bit_parser] requested integer does not fit in the return type");
    if (available < l) fill_buf();
    
    if (l != 64) {
        val = buffer & ((uint64_t(1) << l) - 1);
        buffer >>= l;
    } else {
        val = buffer;
    }
    available -= l;
    idx += l;
    return val;
}

/* skip all zeros from the current position and return the number of skipped zeros */
std::size_t 
ext_bit_parser::parse_0()
{
    std::size_t zeros = 0;
    while (buffer == 0) {
        idx += available;
        zeros += available;
        fill_buf();
    }
    auto l = bit::lsbll(buffer);
    buffer >>= l;
    buffer >>= 1;
    available -= l + 1;
    idx += l + 1;
    return zeros + l;
}

std::size_t 
ext_bit_parser::next_1()
{
    std::optional<std::size_t> pos_in_word = std::nullopt;
    auto buf = buffer;
    while (not (pos_in_word = bit::lsb(buf))) {
        idx += ::bit::size(buf);
        buf = _data[idx >> bit::lsbll(::bit::size(buf))]; // lsbll(64) should be 6
    }
    buffer = buf & (buf - 1);  // clear LSB
    idx = (idx & ~static_cast<std::size_t>(::bit::size(buf) - 1)) + *pos_in_word;
    return idx;
}

std::size_t
ext_bit_parser::get_bit_index() const noexcept
{
    return idx;
}

void 
ext_bit_parser::reset(std::size_t nidx)
{
    if (nidx >= 64 * block_size) throw std::out_of_range("[bit::ext_bit_parser] index out of range");
    idx = nidx;
    buffer = 0;
    available = 0;
}

void 
ext_bit_parser::reset_and_clear_low_bits(std::size_t nidx)
{
    if (nidx >= 64 * block_size) throw std::out_of_range("[bit::ext_bit_parser] index out of range");
    idx = nidx;
    buffer = _data[idx / 64];
    buffer &= uint64_t(-1) << (idx & (64 - 1));  // clear low bits
}

void 
ext_bit_parser::fill_buf() 
{
    buffer = get_next_block();
    available = 64;
}

uint64_t 
ext_bit_parser::get_next_block() const
{
    std::size_t block = idx / 64;
    std::size_t shift = idx % 64;
    if (idx >= 64 * block_size) throw std::out_of_range("[bit::ext_bit_parser] index out of range");
    uint64_t word = _data[block] >> shift;
    if (shift and block + 1 < block_size) word |= _data[block + 1] << (64 - shift);
    return word;
}