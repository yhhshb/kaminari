#include "../include/query/colorsets_accessor.hpp"

colorsets_accessor::colorsets_accessor(){};

colorsets_accessor::colorsets_accessor(colorsets* parent)
    {
        m_parent = parent;
        m_parser = ext_bit_parser(parent->m_mapped_file, 0);
    }

void colorsets_accessor::move_to(uint64_t cid){
    uint64_t bit_index = m_parent->m_offsets.at(cid);
    m_begin = bit_index;
    m_orig = bit_index;
    m_parser.reset(m_begin);

    m_pos_in_list = 0;
    m_prev_val = -1; //2^32 - 1
    m_curr_val = 0;

    m_comp_list_size = 0;
    m_pos_in_comp_list = 0;
    m_comp_val = -1; //2^32 - 1

    m_size = delta_decoder(m_parser);
    if (m_size < m_parent->m_sparse_set_threshold_size) {
        m_type = list_type::delta_gaps;
        m_curr_val = delta_decoder(m_parser);
    } else if (m_size < m_parent->m_dense_set_threshold_size) {
        m_type = list_type::bitmap;
        m_begin = m_parser.get_bit_index();  // after m_size
        m_parser.reset_and_clear_low_bits(m_begin);
        uint64_t pos = m_parser.next_1();
        assert(pos >= m_begin);
        m_curr_val = pos - m_begin;
    } else {
        m_type = list_type::complementary_delta_gaps;
        m_comp_list_size = m_parent->m_num_docs - m_size;
        if (m_comp_list_size > 0) m_comp_val = delta_decoder(m_parser);
        find_next();
    }
}

uint32_t colorsets_accessor::value() const 
{
    return m_curr_val;
}

void colorsets_accessor::next() 
{
    if (m_type == list_type::complementary_delta_gaps) {
        ++m_curr_val;
        if (m_curr_val >= m_parent->m_num_docs) {  // saturate
            m_curr_val = m_parent->m_num_docs;
            return;
        }
        find_next();
    } else if (m_type == list_type::delta_gaps) {
        m_pos_in_list += 1;
        if (m_pos_in_list >= m_size) {  // saturate
            m_curr_val = m_parent->m_num_docs;
            return;
        }
        m_prev_val = m_curr_val;
        m_curr_val = delta_decoder(m_parser) + (m_prev_val + 1);
    } else {
        m_pos_in_list += 1;
        if (m_pos_in_list >= m_size) {  // saturate
            m_curr_val = m_parent->m_num_docs;
            return;
        }
        uint64_t pos = m_parser.next_1();
        assert(pos >= m_begin);
        m_curr_val = pos - m_begin;
    }
}

/* update the state of the iterator to the element which is greater-than or equal-to lower_bound */
void colorsets_accessor::next_geq(uint64_t lower_bound) 
{
    assert(lower_bound <= m_parent->m_num_docs);
    if (m_type == list_type::complementary_delta_gaps) {
        comp_next_geq(lower_bound);
        if (m_curr_val < lower_bound){
            m_curr_val = lower_bound;
        }
        find_next();
    } else {
        while (value() < lower_bound) next();
    }
}

void colorsets_accessor::comp_next_geq(uint64_t lower_bound) 
{
    while (m_comp_val < lower_bound) {
        ++m_pos_in_comp_list;
        if (m_pos_in_comp_list >= m_comp_list_size) break;
        m_prev_val = m_comp_val;
        m_comp_val = delta_decoder(m_parser) + (m_prev_val + 1);
    }
}


uint64_t colorsets_accessor::size() const 
{ 
    return m_size; 
}

colorsets_accessor::list_type colorsets_accessor::type() const 
{ 
    return m_type;
}

uint64_t colorsets_accessor::limit() 
{
    return m_parent->m_num_docs;
}

void colorsets_accessor::find_next() 
{
    while (m_curr_val == m_comp_val) {
        ++m_curr_val;
        ++m_pos_in_comp_list;
        if (m_pos_in_comp_list >= m_comp_list_size) break;
        m_prev_val = m_comp_val;
        m_comp_val = delta_decoder(m_parser) + (m_prev_val + 1);
    }
}