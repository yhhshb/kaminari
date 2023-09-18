#include <string>
// #include <iostream>
#include <numeric>
#include <cassert>

#include "../include/hybrid.hpp"

namespace kaminari {
namespace color_classes {

hybrid::hybrid() 
    : m_num_docs(0), 
      m_sparse_set_threshold_size(0), 
      m_very_dense_set_threshold_size(0),
      m_offsets(dummy_itr_t(0), 0 , 0)
{}

hybrid::builder::builder(std::size_t number_of_documents, bool verb)
    : m_num_docs(number_of_documents), m_num_lists(0), m_num_total_integers(0), verbose(verb)
{
    // if list contains < sparse_set_threshold_size ints, code it with gaps+delta
    m_sparse_set_threshold_size = 0.25 * m_num_docs;

    // if list contains > very_dense_set_threshold_size ints, code it as a complementary set with gaps+delta
    m_very_dense_set_threshold_size = 0.75 * m_num_docs;

    // otherwise: code it as a bitmap of m_num_docs bits
    // m_bvb.reserve(8 * essentials::GB);
    m_offsets.push_back(0);

    if (verbose) {
        std::cerr << "m_num_docs: " << m_num_docs << "\n"
                << "m_sparse_set_threshold_size " << m_sparse_set_threshold_size << "\n" 
                << "m_very_dense_set_threshold_size " << m_very_dense_set_threshold_size << "\n";
    }
}

void 
hybrid::builder::add_color_set(uint32_t* const colors, uint64_t list_size) 
{
    /* encode list_size */
    std::cerr << "before delta (0)\n";
    bit::encoder::delta(m_bvb, list_size);
    std::cerr << "after delta (0)\n";
    if (list_size < m_sparse_set_threshold_size) {
        uint32_t prev_val = colors[0];
        std::cerr << "before delta (1)\n";
        bit::encoder::delta(m_bvb, prev_val);
        std::cerr << "before delta (1)\n";
        for (std::size_t i = 1; i != list_size; ++i) {
            uint32_t val = colors[i];
            assert(val >= prev_val + 1);
            std::cerr << "before delta (2)\n";
            bit::encoder::delta(m_bvb, val - (prev_val + 1));
            std::cerr << "before delta (2)\n";
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
            std::cerr << "before set\n";
            m_bvb.set(old_size + colors[i]);
            std::cerr << "after set\n";
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
                    std::cerr << "before delta (3)\n";
                    bit::encoder::delta(m_bvb, val);
                    std::cerr << "before delta (3)\n";
                    first = false;
                    ++written;
                } else {
                    assert(val >= prev_val + 1);
                    std::cerr << "before delta (4)\n";
                    bit::encoder::delta(m_bvb, val - (prev_val + 1));
                    std::cerr << "before delta (4)\n";
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
            std::cerr << "before delta (5)\n";
            bit::encoder::delta(m_bvb, val - (prev_val + 1));
            std::cerr << "before delta (5)\n";
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
    if (verbose and (m_num_lists % 500000 == 0)) std::cerr << "  processed " << m_num_lists << " lists" << std::endl;
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

    if (verbose) {
        std::cerr << "processed " << m_num_lists << " lists\n" 
                << "m_num_total_integers " << m_num_total_integers << "\n" 
                << "  total bits for ints = " << h.m_colors.size() << "\n" 
                << "  total bits per offsets = " << h.m_offsets.bit_size() << "\n" 
                << "  total bits = " << h.m_offsets.bit_size() + h.m_colors.size() << "n" 
                << "  offsets: " << static_cast<double>(h.m_offsets.bit_size()) / m_num_total_integers << " bits/int\n" 
                << "  lists: " << static_cast<double>(h.m_colors.size()) / m_num_total_integers << " bits/int\n";
    }
}

hybrid::iterator::iterator(hybrid const& ptr, uint64_t begin)
    : m_ptr(ptr)
    , m_begin(begin)
    , m_num_docs(ptr.m_num_docs)
    , m_pos_in_list(0)
    , m_pos_in_comp_list(0)
    , m_comp_list_size(0)
    , m_comp_val(-1)
    , m_prev_val(-1)
    , m_curr_val(0) 
{
    m_parser = bit_parser(ptr.m_colors.data(), ptr.m_colors.block_size(), begin);
    m_size = bit::decoder::delta(m_parser);
    /* set m_type and read the first value */
    if (m_size < ptr.m_sparse_set_threshold_size) {
        m_type = list_type::delta_gaps;
        m_curr_val = bit::decoder::delta(m_parser);
    } else if (m_size < ptr.m_very_dense_set_threshold_size) {
        m_type = list_type::bitmap;
        m_begin = m_parser.get_bit_index();  // after m_size
        m_parser.reset_and_clear_low_bits(m_begin);
        uint64_t pos = m_parser.next_1();
        assert(pos >= m_begin);
        m_curr_val = pos - m_begin;
    } else {
        m_type = list_type::complementary_delta_gaps;
        m_comp_list_size = m_num_docs - m_size;
        if (m_comp_list_size > 0) m_comp_val = bit::decoder::delta(m_parser);
        next_comp_val();
    }
}

// this is needed to annul the next_comp_val() done in the constructor if we want to iterate through the complemented set
void 
hybrid::iterator::reinit_for_complemented_set_iteration() 
{
    assert(m_type == list_type::complementary_delta_gaps);
    m_pos_in_comp_list = 0;
    m_prev_val = -1;
    m_curr_val = 0;
    m_parser = bit_parser(m_ptr.m_colors.data(), m_ptr.m_colors.block_size(), m_begin);
    bit::decoder::delta(m_parser);  // skip m_size
    if (m_comp_list_size > 0) {
        m_comp_val = bit::decoder::delta(m_parser);
    } else {
        m_comp_val = m_num_docs;
    }
}

uint64_t hybrid::iterator::value() const 
{ 
    return m_curr_val; 
}

uint64_t hybrid::iterator::comp_value() const 
{ 
    return m_comp_val; 
}

uint64_t hybrid::iterator::operator*() const 
{ 
    return value(); 
}

const hybrid::iterator&
hybrid::iterator::operator++() 
{
    if (m_type == list_type::complementary_delta_gaps) {
        ++m_curr_val;
        if (m_curr_val >= m_num_docs) {  // saturate
            m_curr_val = m_num_docs;
            return *this;
        }
        next_comp_val();
    } else if (m_type == list_type::delta_gaps) {
        m_pos_in_list += 1;
        if (m_pos_in_list >= m_size) {  // saturate
            m_curr_val = m_num_docs;
            return *this;
        }
        m_prev_val = m_curr_val;
        m_curr_val = bit::decoder::delta(m_parser) + (m_prev_val + 1);
    } else {
        assert(m_type == list_type::bitmap);
        m_pos_in_list += 1;
        if (m_pos_in_list >= m_size) {  // saturate
            m_curr_val = m_num_docs;
            return *this;
        }
        uint64_t pos = m_parser.next_1();
        assert(pos >= m_begin);
        m_curr_val = pos - m_begin;
    }
    return *this;
}

void 
hybrid::iterator::next_comp() 
{
    ++m_pos_in_comp_list;
    if (m_pos_in_comp_list >= m_comp_list_size) {  // saturate
        m_comp_val = m_num_docs;
        return;
    }
    m_prev_val = m_comp_val;
    m_comp_val = bit::decoder::delta(m_parser) + (m_prev_val + 1);
}

/* update the state of the iterator to the element
    which is greater-than or equal-to lower_bound */
void 
hybrid::iterator::next_geq(uint64_t lower_bound) 
{
    assert(lower_bound <= m_num_docs);
    if (m_type == list_type::complementary_delta_gaps) {
        next_geq_comp_val(lower_bound);
        m_curr_val = lower_bound;
    } else {
        while (value() < lower_bound) m_parser.next_1();
    }
    assert(value() >= lower_bound);
}

std::size_t 
hybrid::iterator::size() const 
{ 
    return m_size; 
}

std::size_t 
hybrid::iterator::num_docs() const 
{ 
    return m_num_docs; 
}

hybrid::list_type 
hybrid::iterator::type() const 
{ 
    return m_type;
}

void 
hybrid::iterator::next_comp_val() 
{
    while (m_curr_val == m_comp_val) {
        ++m_curr_val;
        ++m_pos_in_comp_list;
        if (m_pos_in_comp_list >= m_comp_list_size) break;
        m_prev_val = m_comp_val;
        m_comp_val = bit::decoder::delta(m_parser) + (m_prev_val + 1);
    }
}

void 
hybrid::iterator::next_geq_comp_val(uint64_t lower_bound) 
{
    while (m_comp_val < lower_bound) {
        ++m_pos_in_comp_list;
        if (m_pos_in_comp_list >= m_comp_list_size) break;
        m_prev_val = m_comp_val;
        m_comp_val = bit::decoder::delta(m_parser) + (m_prev_val + 1);
    }
}

hybrid::iterator 
hybrid::colors(uint64_t color_class_id) const 
{
    assert(color_class_id < num_color_classes());
    // uint64_t begin = m_offsets.access(color_class_id);
    uint64_t begin = m_offsets.at(color_class_id);
    return iterator(*this, begin);
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
        m_colors.size() + //essentials::vec_bytes(m_colors) +
        m_offsets.bit_size();
}

void 
hybrid::print_stats() const 
{
    uint64_t num_buckets = 10;
    assert(num_buckets > 0);
    uint64_t bucket_size = m_num_docs / num_buckets;
    std::vector<uint32_t> list_size_upperbounds;
    for (std::size_t i = 0, curr_list_size_upper_bound = bucket_size; 
        i != num_buckets; 
        ++i, curr_list_size_upper_bound += bucket_size) 
    {
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
        // uint64_t offset = m_offsets.access(color_class_id);
        uint64_t offset = m_offsets.at(color_class_id);
        bit_parser it(m_colors.data(), m_colors.block_size(), offset);
        uint32_t list_size = bit::decoder::delta(it);
        // uint64_t num_bits = m_offsets.access(color_class_id + 1) - offset;
        uint64_t num_bits = m_offsets.at(color_class_id + 1) - offset;
        auto bucket_it = std::upper_bound(list_size_upperbounds.begin(), list_size_upperbounds.end(), list_size);
        if (bucket_it != list_size_upperbounds.begin() and *(bucket_it - 1) == list_size) --bucket_it;
        uint64_t bucket_index = std::distance(list_size_upperbounds.begin(), bucket_it);
        num_bits_per_bucket[bucket_index] += num_bits;
        num_lists_per_bucket[bucket_index] += 1;
        num_ints_per_bucket[bucket_index] += list_size;
        num_total_integers += list_size;
    }

    if (verbose) std::cerr << "CCs SPACE BREAKDOWN:\n";
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
            if (verbose) {
                std::cerr   << "num. lists of size > " 
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
    }
    assert(integers == num_total_integers);
    assert(std::accumulate(num_lists_per_bucket.begin(), num_lists_per_bucket.end(), uint64_t(0)) == num_lists);
    if (verbose) {
        std::cerr << "  colors: " << static_cast<double>(bits) / integers << " bits/int\n"
                << "  offsets: "
                << static_cast<double>((sizeof(m_num_docs) + sizeof(m_sparse_set_threshold_size) + sizeof(m_very_dense_set_threshold_size)) * 8 + m_offsets.bit_size()) / integers
                << " bits/int\n";
    }
}

} // namespace color_classes 
} // namespace kaminari