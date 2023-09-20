#ifndef KAMINARI_MAPPER_ARRAY_HPP
#define KAMINARI_MAPPER_ARRAY_HPP

#include <map>

#include "../bundled/biolib/include/kmer_view.hpp"
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
        uint32_t at(std::size_t kmer_idx) const;
        std::map<std::size_t, std::size_t> get_histogram() const noexcept;
        const std::vector<std::size_t>& vector_data() const noexcept {return map;}

        template <class Visitor>
        void visit(Visitor& visitor) const;

        template <class Visitor>
        void visit(Visitor& visitor);

    private:
        std::vector<std::size_t> map;
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
            // try {
                if (!same_color) {
                    ++color_class_id; // color_id
                }
                auto hash_values = hf(unitig.data, unitig.size, true);
                for (auto v : hash_values) uncompressed_map[v] = color_class_id;
            // } catch (std::exception const& e) {
            //     std::cerr << e.what() << std::endl;
            //     exit(1);
            // }
        }
    );
    
    map.swap(uncompressed_map); // FIXME with proper compression
}

template <class MPHF>
uint32_t
array_based<MPHF>::at(std::size_t kmer_idx) const
{
    return map.at(kmer_idx);
}

template <class MPHF>
std::map<std::size_t, std::size_t> 
array_based<MPHF>::get_histogram() const noexcept
{
    std::map<std::size_t, std::size_t> hist;
    for (auto v : map) {
        if (hist.find(v) == hist.end()) hist[v] = 1;
        else hist[v] += 1;
    } 
    return hist;
}

template <class MPHF>
template <class Visitor>
void 
array_based<MPHF>::visit(Visitor& visitor) const
{
    visitor.visit(map);
}

template <class MPHF>
template <class Visitor>
void 
array_based<MPHF>::visit(Visitor& visitor)
{
    visitor.visit(map);
}

} // namespace mapper
} // namespace kaminari

#endif // KAMINARI_MAPPER_ARRAY_HPP
