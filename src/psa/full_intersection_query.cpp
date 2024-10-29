#include "../../include/index.hpp"
#include "../../include/hybrid.hpp"

namespace kaminari {
namespace minimizer {

#define CLASS_HEADER template <class ColorClasses, class ColorMapper>
#define METHOD_HEADER index<ColorClasses, ColorMapper>

template class index<kaminari::color_classes::hybrid, pthash::compact_vector>;

CLASS_HEADER
std::vector<typename METHOD_HEADER::color_t> 
METHOD_HEADER::query_full_intersection(char const * const q, const std::size_t l, options_t& opts) const noexcept
{
    if (opts.verbose > 2) std::cerr << "step 1: collect color class ids\n";
    std::vector<std::size_t> ccids;
    { // collect color class ids
        std::size_t contig_mmer_count;
        std::vector<::minimizer::record_t> mms_buffer;
        [[maybe_unused]] auto contig_kmer_count = ::minimizer::from_string<hash64>(q, l, k, m, seed, canonical, contig_mmer_count, mms_buffer);
        for (const auto& record : mms_buffer) { 
            ccids.push_back(m_map[(hf(record.itself))]);
        }
    }

    if (opts.verbose > 2) std::cerr << "step 2: ids to colors\n";
    std::vector<typename ColorClasses::row_accessor> color_itrs;
    bool all_very_dense = true;
    {
        auto last = std::unique(ccids.begin(), ccids.end()); // deduplicate color class ids
        for (auto itr = ccids.begin(); itr != last; ++itr) { 
            color_itrs.push_back(m_ccs.colors_at(*itr));
            if (color_itrs.back().type() != ColorClasses::row_accessor::complementary_delta_gaps) {
                all_very_dense = false;
            }
        }
    }
    if (opts.verbose > 2) {
        std::cerr << "step 3: computing intersections\n";
        if (all_very_dense) std::cerr << "\tcompute dense intersection\n";
        else std::cerr << "\tcompute mixed intersection (dense and sparse vectors)\n"; 
    }
    if (color_itrs.empty()) return {};
    if (all_very_dense) return full_dense_intersection(std::move(color_itrs)); // intersect of dense rows
    else return full_mixed_intersection(std::move(color_itrs)); // intersect dense and sparse rows
}


CLASS_HEADER
std::vector<typename METHOD_HEADER::color_t> 
METHOD_HEADER::full_dense_intersection(std::vector<typename ColorClasses::row_accessor>&& color_id_itrs) const noexcept
{
    std::vector<color_t> colors;
    std::size_t vec_size = color_id_itrs.size();
    std::size_t filenames_size = m_filenames.size();

    std::vector<bool> presence(filenames_size, true);

    for (uint64_t i = 0; i != vec_size; ++i) {
        while (color_id_itrs[i].comp_value() != filenames_size) {
            presence[color_id_itrs[i].comp_value()] = false;
            color_id_itrs[i].comp_next();
        }
    }

    for (uint32_t i = 0; i != filenames_size; ++i) {
        if (presence[i]) colors.push_back(i);
    }

    return colors;
}

CLASS_HEADER
std::vector<typename METHOD_HEADER::color_t> 
METHOD_HEADER::full_mixed_intersection(std::vector<typename ColorClasses::row_accessor>&& color_id_itrs) const noexcept
{
    std::vector<color_t> colors;
    std::size_t vec_size = color_id_itrs.size();
    std::size_t filenames_size = m_filenames.size();

    std::vector<uint32_t> counts(filenames_size, 0);

    for (uint64_t i = 0; i != vec_size; ++i) {
        while (color_id_itrs[i].value() != filenames_size) {
            counts[color_id_itrs[i].value()] += 1;
            color_id_itrs[i].next();
        }
    }

    for (uint32_t i = 0; i != filenames_size; ++i) {
        if (counts[i] == vec_size) colors.push_back(i);
    }

    return colors;
}


#undef CLASS_HEADER
#undef METHOD_HEADER

} // namespace minimizer
} // namespace kaminari