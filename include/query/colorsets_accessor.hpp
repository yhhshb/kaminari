#ifndef COLORSETS_ACCESSOR_HPP
#define COLORSETS_ACCESSOR_HPP

#include "../codes.hpp"
#include "../mmap.hpp"
#include "../ext_bit_parser.hpp"
#include "../../bundled/biolib/include/elias_fano.hpp"
#include "../../bundled/biolib/include/io.hpp"

class colorsets_accessor
    // class used to access the colors contained in a color class, works like a reader head, can be re-located elsewhere to read another colorset
        {
            public:
                enum list_type { 
                    delta_gaps = 0,
                    bitmap = 1,
                    complementary_delta_gaps = 2
                };

                colorsets_accessor();
                colorsets_accessor(std::string basename);
                void move_to(uint64_t cid);
                uint32_t value() const;
                void next();
                void next_geq(uint64_t lower_bound);
                void comp_next_geq(uint64_t lower_bound);

                void reset();

                list_type type() const;
                uint64_t size() const;
                uint64_t limit();
                
            private:
                uint64_t m_num_docs;
                uint64_t m_sparse_set_threshold_size;
                uint64_t m_dense_set_threshold_size;
                bit::ef::array m_offsets;

                ext_bit_parser m_parser;
                uint64_t m_begin;
                uint64_t m_orig;
                list_type m_type;

                uint32_t m_size;
                uint32_t m_pos_in_list;
                uint32_t m_prev_val;
                uint32_t m_curr_val;
                
                uint32_t m_comp_list_size;
                uint32_t m_pos_in_comp_list;
                uint32_t m_comp_val;
                
                void find_next();
        };

#endif // COLORSETS_ACCESSOR_HPP