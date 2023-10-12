#ifndef KAMINARI_MAPPER_ARRAY_HPP
#define KAMINARI_MAPPER_ARRAY_HPP

#include <map>
#include <cassert>
#include "../bundled/biolib/include/rle_view.hpp"
#include "../bundled/biolib/include/packed_vector.hpp"
#include "../bundled/biolib/include/elias_fano.hpp"
#include "../bundled/biolib/include/cumulative_iterator.hpp"
#include "constants.hpp"
#include "GGCAT.hpp"

namespace kaminari {
namespace mapper {

template <class MPHF>
class array_based 
{
    public:
        array_based() {}
        void build(MPHF const& hf, GGCAT& cdbg);
        uint32_t at(std::size_t kmer_idx) const {return map.at(kmer_idx);}
        // std::map<std::size_t, std::size_t> get_histogram() const noexcept;

        template <class Visitor> void visit(Visitor& visitor) const {visitor.visit(map);}
        template <class Visitor> void visit(Visitor& visitor) {visitor.visit(map);}

    private:
        class rle_encoding
        {
            public:
                using bv_t = bit::packed::vector<std::size_t>;
                using ef_t = bit::ef::array;
                rle_encoding() noexcept : values(bv_t(0)) {}
                void build(std::vector<std::size_t> const& uncompressed_map);
                uint32_t at(std::size_t idx) const;
                std::size_t size() const noexcept {return values.size();}
                template <class Visitor> void visit(Visitor& visitor) const {
                    visitor.visit(values);
                    visitor.visit(lengths);
                }
                template <class Visitor> void visit(Visitor& visitor) {
                    visitor.visit(values);
                    visitor.visit(lengths);
                }
            private:
                bv_t values;
                ef_t lengths;
        };
        // std::vector<std::size_t> 
        rle_encoding map;
};

template <class MPHF>
void 
array_based<MPHF>::build(MPHF const& hf, GGCAT& cdbg) 
{
    std::size_t color_class_id = 0;
    std::vector<std::size_t> uncompressed_map;
    uncompressed_map.resize(hf.get_kmer_count());
    
    cdbg.loop_through_unitigs(
        [&](ggcat::Slice<char> const unitig, ggcat::Slice<uint32_t> const colors, bool same_color) 
        {
            if (!same_color) {
                ++color_class_id; // color_id
            }
            auto hash_values = hf(unitig.data, unitig.size, true);
            for (auto v : hash_values) uncompressed_map[v] = color_class_id;
        }
    );
    // map.swap(uncompressed_map); // FIXME with proper compression
    map.build(uncompressed_map);
}

// template <class MPHF>
// std::map<std::size_t, std::size_t> 
// array_based<MPHF>::get_histogram() const noexcept
// {
//     std::map<std::size_t, std::size_t> hist;
//     for (auto v : map) {
//         if (hist.find(v) == hist.end()) hist[v] = 1;
//         else hist[v] += 1;
//     } 
//     return hist;
// }

template <class MPHF>
void
array_based<MPHF>::rle_encoding::build(std::vector<std::size_t> const& uncompressed_map)
{
    wrapper::rle_view view(uncompressed_map.cbegin(), uncompressed_map.cend());
    std::vector<std::size_t> uncompressed_lengths;
    std::size_t length_check = 0;
    for (auto itr = view.cbegin(); itr != view.cend(); ++itr) {
        auto [val, len] = *itr;
        values.push_back(val);
        uncompressed_lengths.push_back(len); // IDEA: re-use uncompressed_map to avoid another vector
        length_check += len;
    }
    if (length_check != uncompressed_map.size()) throw std::logic_error("Manually calculated number of values != vector size");
    lengths = ef_t(iterators::cumulative_iterator(uncompressed_lengths.cbegin()), uncompressed_lengths.size(), uncompressed_map.size());
    assert(values.size() == lengths.size());
}

template <class MPHF>
uint32_t
array_based<MPHF>::rle_encoding::at(std::size_t idx) const
{
    std::size_t i = lengths.leq_find(idx); // find the index the right-most smallest than idx cumulative length
    return values.template at<uint32_t>(i);
}

} // namespace mapper
} // namespace kaminari

#endif // KAMINARI_MAPPER_ARRAY_HPP
