#ifndef EXT_BIT_PARSER_HPP
#define EXT_BIT_PARSER_HPP

#include "../bundled/biolib/include/bit_operations.hpp"
#include "mmap.hpp"

class ext_bit_parser
{
    public:
        ext_bit_parser();
        ext_bit_parser(uint64_t const * const data, std::size_t size, std::size_t starting_bit_position = 0);
        ext_bit_parser(mymm::immap<uint64_t>& mapped_file, std::size_t starting_bit_position = 0);

        uint64_t parse_fixed(std::size_t l);
        std::size_t parse_0();
        std::size_t next_1();

        std::size_t get_bit_index() const noexcept;
        void reset(std::size_t nidx);
        void reset_and_clear_low_bits(std::size_t nidx);
        
    private:
        void fill_buf();
        uint64_t get_next_block() const;

        uint64_t const* _data;
        std::size_t block_size;
        std::size_t idx;
        uint64_t buffer;
        std::size_t available; // buffer handling
};

#endif // EXT_BIT_PARSER_HPP