#ifndef COLORSETS_ACCESSOR_HPP
#define COLORSETS_ACCESSOR_HPP

#include "../codes.hpp"
#include "../mmap.hpp"
#include "../ext_bit_parser.hpp"
#include "../colorsets.hpp"
#include "../../bundled/biolib/include/elias_fano.hpp"
#include "../../bundled/biolib/include/io.hpp"

class colorsets_accessor
    // class used to access the colors contained in a color class, works like a reader head, can be re-located elsewhere to read another colorset
    // idea: 1 per query thread
        {
            public:
                enum list_type { 
                    delta_gaps = 0,
                    bitmap = 1,
                    complementary_delta_gaps = 2
                };

                colorsets_accessor();
                colorsets_accessor(colorsets* parent);
                void move_to(uint64_t cid);
                uint32_t value() const;
                void next();
                void next_geq(uint64_t lower_bound);
                void comp_next_geq(uint64_t lower_bound);

                list_type type() const;
                uint64_t size() const;
                uint64_t limit();
                
            private:
                colorsets* m_parent;   // non owning so that each thread uses the same offsets elias fano
                ext_bit_parser m_parser;     // lightweight parser into mmap
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