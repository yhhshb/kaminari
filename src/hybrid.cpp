#include <string>
#include <numeric>
#include <iostream>
#include <cassert>

#include "../include/hybrid.hpp"
#include "../bundled/biolib/include/counting_iterator.hpp"

// #include <iostream>

namespace kaminari {
namespace color_classes {

typedef iterators::counting_iterator<uint32_t> dummy_itr_t;

hybrid::hybrid() 
    : m_num_docs(0), 
      m_sparse_set_threshold_size(0), 
      m_very_dense_set_threshold_size(0),
      m_offsets(dummy_itr_t(0), 0 , 0)
{}

hybrid::builder::builder(std::size_t number_of_documents, std::size_t verbosity_level)
    : m_num_docs(number_of_documents), m_num_lists(0), m_num_total_integers(0), m_verbosity_level(verbosity_level)
{
    // if list contains < sparse_set_threshold_size ints, code it with gaps+delta
    m_sparse_set_threshold_size = 0.25 * m_num_docs;

    // if list contains > very_dense_set_threshold_size ints, code it as a complementary set with gaps+delta
    m_very_dense_set_threshold_size = 0.75 * m_num_docs;

    // otherwise: code it as a bitmap of m_num_docs bits
    // m_bvb.reserve(8 * essentials::GB);
    m_offsets.push_back(0);

    if (verbosity_level > 0) {
        std::cerr 
            << "m_num_docs: " << m_num_docs << "\n"
            << "m_sparse_set_threshold_size " << m_sparse_set_threshold_size << "\n" 
            << "m_very_dense_set_threshold_size " << m_very_dense_set_threshold_size << "\n";
    }
}

void 
hybrid::builder::add_color_set(color_t const * const colors, std::size_t list_size) 
{
    /* encode list_size */
    bit::encoder::delta(m_bvb, list_size);
    if (list_size < m_sparse_set_threshold_size) {
        uint32_t prev_val = colors[0];
        bit::encoder::delta(m_bvb, prev_val);
        for (std::size_t i = 1; i != list_size; ++i) {
            uint32_t val = colors[i];
            assert(val >= prev_val + 1);
            bit::encoder::delta(m_bvb, val - (prev_val + 1));
            prev_val = val;
        }
    } else if (list_size < m_very_dense_set_threshold_size) {
        // bit_vector bvb_ints;
        // bvb_ints.resize(m_num_docs);
        // for (std::size_t i = 0; i != list_size; ++i) bvb_ints.set(colors[i]);
        // m_bvb.append(bvb_ints);
        auto old_size = m_bvb.size();
        m_bvb.resize(old_size + m_num_docs);
        for (std::size_t i = 0; i != list_size; ++i) {
            m_bvb.set(old_size + colors[i]);
        }
    } else {
        bool first = true;
        uint32_t val = 0;
        uint32_t prev_val = -1;
        uint32_t written = 0;
        for (std::size_t i = 0; i != list_size; ++i) {
            uint32_t x = colors[i];
            while (val < x) {
                if (first) {
                    bit::encoder::delta(m_bvb, val);
                    first = false;
                    ++written;
                } else {
                    assert(val >= prev_val + 1);
                    bit::encoder::delta(m_bvb, val - (prev_val + 1));
                    ++written;
                }
                prev_val = val;
                ++val;
            }
            assert(val == x);
            val++;  // skip x
        }
        while (val < m_num_docs) {
            assert(val >= prev_val + 1);
            bit::encoder::delta(m_bvb, val - (prev_val + 1));
            prev_val = val;
            ++val;
            ++written;
        }
        assert(val == m_num_docs);
        /* complementary_list_size = m_num_docs - list_size */
        assert(m_num_docs - list_size <= m_num_docs);
        assert(written == m_num_docs - list_size);
    }
    m_offsets.push_back(m_bvb.size());
    m_num_total_integers += list_size;
    m_num_lists += 1;
    if (m_verbosity_level > 0 and (m_num_lists % 500000 == 0)) {
        std::cerr << "processed " << m_num_lists << " lists\n";
    }
}

void 
hybrid::builder::build(hybrid& h) 
{
    h.m_num_docs = m_num_docs;
    h.m_sparse_set_threshold_size = m_sparse_set_threshold_size;
    h.m_very_dense_set_threshold_size = m_very_dense_set_threshold_size;
    assert(m_num_lists == m_offsets.size() - 1);

    // h.m_offsets.encode(m_offsets.begin(), m_offsets.size(), m_offsets.back());
    h.m_offsets = ef_sequence(m_offsets.begin(), m_offsets.size(), m_offsets.back());
    h.m_colors.swap(m_bvb);

    if (m_verbosity_level > 0) {
        std::cerr 
            << "processed " << m_num_lists << " lists\n" 
            << "\tm_num_total_integers " << m_num_total_integers << "\n" 
            << "\ttotal bits for ints = " << h.m_colors.size() << "\n" 
            << "\ttotal bits per offsets = " << h.m_offsets.bit_size() << "\n" 
            << "\ttotal bits = " << h.m_offsets.bit_size() + h.m_colors.size() << "n" 
            << "\toffsets: " << static_cast<double>(h.m_offsets.bit_size()) / m_num_total_integers << " bits/int\n" 
            << "\tlists: " << static_cast<double>(h.m_colors.size()) / m_num_total_integers << " bits/int\n";
    }
}

hybrid::row_accessor::row_accessor(hybrid const* parent_colors_storage, std::size_t start_idx)
    : m_parent(parent_colors_storage)
    , m_begin(start_idx)
    // , m_num_docs(ptr.m_num_docs)
    , m_pos_in_list(0)
    , m_prev_val(-1)
    , m_curr_val(0) 
    , m_comp_list_size(0)
    , m_pos_in_comp_list(0)
    , m_comp_val(-1)
    
{
    m_parser = bit_parser(m_parent->m_colors.data(), m_parent->m_colors.block_size(), m_begin);
    m_size = bit::decoder::delta(m_parser);
    /* set m_type and read the first value */
    if (m_size < m_parent->m_sparse_set_threshold_size) {
        m_type = list_type::delta_gaps;
        m_curr_val = bit::decoder::delta(m_parser);
    } else if (m_size < m_parent->m_very_dense_set_threshold_size) {
        m_type = list_type::bitmap;
        m_begin = m_parser.get_bit_index();  // after m_size
        m_parser.reset_and_clear_low_bits(m_begin);
        uint64_t pos = m_parser.next_1();
        assert(pos >= m_begin);
        m_curr_val = pos - m_begin;
    } else {
        m_type = list_type::complementary_delta_gaps;
        m_comp_list_size = m_parent->m_num_docs - m_size;
        if (m_comp_list_size > 0) m_comp_val = bit::decoder::delta(m_parser);
        find_next();
    }
}

hybrid::row_accessor::value_type
hybrid::row_accessor::value() const 
{
    return m_curr_val;
}

void
hybrid::row_accessor::next() 
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
        m_curr_val = bit::decoder::delta(m_parser) + (m_prev_val + 1);
    } else {
        assert(m_type == list_type::bitmap);
        m_pos_in_list += 1;
        if (m_pos_in_list >= m_size) {  // saturate
            m_curr_val = m_parent->m_num_docs;
            return;
        }
        uint64_t pos = m_parser.next_1();
        assert(pos >= m_begin);
        m_curr_val = pos - m_begin;
    }
    return;
}

/* update the state of the iterator to the element which is greater-than or equal-to lower_bound */
void 
hybrid::row_accessor::next_geq(uint64_t lower_bound) 
{
    assert(lower_bound <= m_parent->m_num_docs);
    if (m_type == list_type::complementary_delta_gaps) {
        comp_next_geq(lower_bound);
        m_curr_val = lower_bound;
    } else {
        while (value() < lower_bound) m_parser.next_1();
    }
    assert(value() >= lower_bound);
}

// this is needed in order to void the effects of find_next() in the constructor
void 
hybrid::row_accessor::reinit_for_complemented_set_iteration() 
{
    assert(m_type == list_type::complementary_delta_gaps);
    if (m_type != list_type::complementary_delta_gaps) 
        throw std::runtime_error("[row_accessor] only accessors storing complementary_delta_gaps can be transformed");
    m_pos_in_comp_list = 0;
    m_prev_val = -1;
    m_curr_val = 0;
    m_parser = bit_parser(m_parent->m_colors.data(), m_parent->m_colors.block_size(), m_begin);
    bit::decoder::delta(m_parser);  // skip m_size
    if (m_comp_list_size > 0) {
        m_comp_val = bit::decoder::delta(m_parser);
    } else {
        m_comp_val = m_parent->m_num_docs;
    }
}

hybrid::row_accessor::value_type 
hybrid::row_accessor::comp_value() const 
{ 
    return m_comp_val; 
}

void 
hybrid::row_accessor::comp_next() 
{
    ++m_pos_in_comp_list;
    if (m_pos_in_comp_list >= m_comp_list_size) {  // saturate
        m_comp_val = m_parent->m_num_docs;
        return;
    }
    m_prev_val = m_comp_val;
    m_comp_val = bit::decoder::delta(m_parser) + (m_prev_val + 1);
}

void 
hybrid::row_accessor::comp_next_geq(uint64_t lower_bound) 
{
    while (m_comp_val < lower_bound) {
        ++m_pos_in_comp_list;
        if (m_pos_in_comp_list >= m_comp_list_size) break;
        m_prev_val = m_comp_val;
        m_comp_val = bit::decoder::delta(m_parser) + (m_prev_val + 1);
    }
}

std::size_t 
hybrid::row_accessor::size() const 
{ 
    return m_size; 
}

hybrid::row_accessor::list_type 
hybrid::row_accessor::type() const 
{ 
    return m_type;
}

void 
hybrid::row_accessor::find_next() 
{
    while (m_curr_val == m_comp_val) {
        ++m_curr_val;
        ++m_pos_in_comp_list;
        if (m_pos_in_comp_list >= m_comp_list_size) break;
        m_prev_val = m_comp_val;
        m_comp_val = bit::decoder::delta(m_parser) + (m_prev_val + 1);
    }
}

// std::size_t 
// hybrid::iterator::num_docs() const 
// { 
//     return m_num_docs; 
// }

hybrid::row_accessor 
hybrid::colors_at(std::size_t color_class_id) const 
{
    assert(color_class_id < num_color_classes());
    uint64_t begin = m_offsets.at(color_class_id);
    return row_accessor(this, begin);
}

std::size_t 
hybrid::num_docs() const 
{ 
    return m_num_docs;
}

std::size_t 
hybrid::num_color_classes() const 
{
    return m_offsets.size() - 1;
}

std::size_t 
hybrid::num_bits() const 
{
    return 
        8 * (sizeof(m_num_docs) + sizeof(m_sparse_set_threshold_size) + sizeof(m_very_dense_set_threshold_size)) + 
        m_colors.size() + 
        m_offsets.bit_size();
}

void 
hybrid::print_stats(std::ostream& out) const 
{
    uint64_t num_buckets = 10;
    assert(num_buckets > 0);
    uint64_t bucket_size = m_num_docs / num_buckets;
    std::vector<uint32_t> list_size_upperbounds;
    for (
        std::size_t i = 0, curr_list_size_upper_bound = bucket_size; 
        i != num_buckets; 
        ++i, curr_list_size_upper_bound += bucket_size
    ) {
        if (i == num_buckets - 1) curr_list_size_upper_bound = m_num_docs;
        list_size_upperbounds.push_back(curr_list_size_upper_bound);
    }

    std::vector<uint64_t> num_bits_per_bucket;
    std::vector<uint64_t> num_lists_per_bucket;
    std::vector<uint64_t> num_ints_per_bucket;
    num_bits_per_bucket.resize(num_buckets, 0);
    num_lists_per_bucket.resize(num_buckets, 0);
    num_ints_per_bucket.resize(num_buckets, 0);

    const uint64_t num_lists = num_color_classes();
    uint64_t num_total_integers = 0;
    for (uint64_t color_class_id = 0; color_class_id != m_offsets.size() - 1; ++color_class_id) {
        uint64_t offset = m_offsets.at(color_class_id);
        bit_parser it(m_colors.data(), m_colors.block_size(), offset);
        uint32_t list_size = bit::decoder::delta(it);
        uint64_t num_bits = m_offsets.at(color_class_id + 1) - offset;
        auto bucket_it = std::upper_bound(list_size_upperbounds.begin(), list_size_upperbounds.end(), list_size);
        if (bucket_it != list_size_upperbounds.begin() and *(bucket_it - 1) == list_size) --bucket_it;
        uint64_t bucket_index = std::distance(list_size_upperbounds.begin(), bucket_it);
        num_bits_per_bucket[bucket_index] += num_bits;
        num_lists_per_bucket[bucket_index] += 1;
        num_ints_per_bucket[bucket_index] += list_size;
        num_total_integers += list_size;
    }

    out << "CCs SPACE BREAKDOWN:\n";
    uint64_t integers = 0;
    uint64_t bits = 0;
    const uint64_t total_bits = num_bits();
    for (std::size_t i = 0, curr_list_size_upper_bound = 0; i != num_buckets; ++i) {
        if (i == num_buckets - 1) curr_list_size_upper_bound = m_num_docs;
        else curr_list_size_upper_bound += bucket_size;
        if (num_lists_per_bucket[i] > 0) {
            uint64_t n = num_ints_per_bucket[i];
            integers += n;
            bits += num_bits_per_bucket[i];
            out   
                << "num. lists of size > " 
                << (curr_list_size_upper_bound - bucket_size)
                << " and <= " 
                << curr_list_size_upper_bound 
                << ": "
                << num_lists_per_bucket[i] 
                << " ("
                << (num_lists_per_bucket[i] * 100.0) / num_lists
                << "%) -- integers: " 
                << n 
                << " (" << (n * 100.0) / num_total_integers
                << "%) -- bits/int: " << static_cast<double>(num_bits_per_bucket[i]) / n
                << " -- "
                << static_cast<double>(num_bits_per_bucket[i]) / total_bits * 100.0
                << "\% of total space" << '\n';
        }
    }
    assert(integers == num_total_integers);
    assert(std::accumulate(num_lists_per_bucket.begin(), num_lists_per_bucket.end(), uint64_t(0)) == num_lists);
    out 
        << "  colors: " << static_cast<double>(bits) / integers << " bits/int\n"
        << "  offsets: "
        << static_cast<double>((sizeof(m_num_docs) + sizeof(m_sparse_set_threshold_size) + sizeof(m_very_dense_set_threshold_size)) * 8 + m_offsets.bit_size()) / integers
        << " bits/int\n";
}

} // namespace color_classes 
} // namespace kaminari