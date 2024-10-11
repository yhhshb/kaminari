#include "../../include/index.hpp"
#include "../../include/hybrid.hpp"

namespace kaminari {
namespace minimizer {

#define CLASS_HEADER template <class ColorClasses, class ColorMapper>
#define METHOD_HEADER index<ColorClasses, ColorMapper>

template class index<kaminari::color_classes::hybrid, pthash::compact_vector>;

CLASS_HEADER
std::vector<typename METHOD_HEADER::color_t> 
METHOD_HEADER::query_full_intersection(char const * const q, const std::size_t l, std::size_t verbosity_level) const noexcept
{
    if (verbosity_level > 0) std::cerr << "step 1: collect color class ids\n";
    std::vector<std::size_t> ccids;
    { // collect color class ids
        std::size_t contig_mmer_count;
        std::vector<::minimizer::record_t> mms_buffer;
        [[maybe_unused]] auto contig_kmer_count = ::minimizer::from_string<hash64>(q, l, k, m, seed, canonical, contig_mmer_count, mms_buffer);
        for (const auto& record : mms_buffer) { 
            ccids.push_back(m_map[(hf(record.itself))]);
        }
    }

    if (verbosity_level > 0) std::cerr << "step 2: ids to colors\n";
    std::vector<typename ColorClasses::row_accessor> color_itrs;
    bool all_very_dense = true;
    {
        auto last = std::unique(ccids.begin(), ccids.end()); // deduplicate color class ids
        for (auto itr = ccids.begin(); itr != last; ++itr) { 
            color_itrs.push_back(m_ccs.colors_at(*itr));

            // IMPROVEMENT 
            // Divide rows based on dense/sparse 
            // 1) apply dense_intersection to the dense dataset 
            // 2) and the complementary iterators from the sparse set <-- Not sure about this
            // make the complementary of the result of (2) and make another intersection for the final result
            if (color_itrs.back().type() != ColorClasses::row_accessor::complementary_delta_gaps) {
                all_very_dense = false;
            }
        }
    }
    if (verbosity_level > 0) {
        std::cerr << "step 3: computing intersections\n";
        if (all_very_dense) std::cerr << "\tcompute dense intersection\n";
        else std::cerr << "\tcompute mixed intersection (dense and sparse vectors)\n"; 
    }
    if (color_itrs.empty()) return {};
    if (all_very_dense) return full_dense_intersection(std::move(color_itrs), verbosity_level); // intersect of dense rows
    else return full_mixed_intersection(std::move(color_itrs), verbosity_level); // intersect dense and sparse rows
}


CLASS_HEADER
std::vector<typename METHOD_HEADER::color_t> 
METHOD_HEADER::full_dense_intersection(std::vector<typename ColorClasses::row_accessor>&& color_id_itrs, std::size_t verbosity_level) const noexcept
{
    std::vector<color_t> tmp;
    { // step 1: take the union of complementary sets
        for (auto& itr : color_id_itrs) itr.reinit_for_complemented_set_iteration();
        color_t candidate = (
            *std::min_element(
                color_id_itrs.begin(), 
                color_id_itrs.end(),
                [](auto const& x, auto const& y) {
                    return x.comp_value() < y.comp_value();
                }
            )
        ).comp_value();
        // const uint32_t num_docs = iterators[0].num_docs();

        tmp.reserve(m_filenames.size());
        while (candidate < m_filenames.size()) {
            color_t next_candidate = m_filenames.size();
            for (uint64_t i = 0; i != color_id_itrs.size(); ++i) {
                if (color_id_itrs.at(i).comp_value() == candidate) color_id_itrs[i].comp_next();
                /* compute next minimum */
                if (color_id_itrs.at(i).comp_value() < next_candidate) {
                    next_candidate = color_id_itrs.at(i).comp_value();
                }
            }
            tmp.push_back(candidate);
            assert(next_candidate > candidate);
            candidate = next_candidate;
        }
    }
    std::vector<color_t> colors;
    { // step 2: compute the intersection by scanning tmp
        color_t candidate = 0;
        for (std::size_t i = 0; i != tmp.size(); ++i) {
            while (candidate < tmp[i]) {
                colors.push_back(candidate);
                candidate += 1;
            }
            candidate += 1;  // skip the candidate because it is equal to tmp[i]
        }
        while (candidate < m_filenames.size()) {
            colors.push_back(candidate);
            candidate += 1;
        }
    }
    return colors;
}

CLASS_HEADER
std::vector<typename METHOD_HEADER::color_t> 
METHOD_HEADER::full_mixed_intersection(std::vector<typename ColorClasses::row_accessor>&& color_id_itrs, std::size_t verbosity_level) const noexcept
{
    std::sort(
        color_id_itrs.begin(),
        color_id_itrs.end(),
        [](auto const& x, auto const& y) { return x.size() < y.size(); }
    );

    std::vector<color_t> colors;
    {
        // const uint32_t num_docs = iterators[0].num_docs();
        color_t candidate = color_id_itrs.front().value();
        std::size_t i = 1;
        while (candidate < m_filenames.size()) {
            for (; i != color_id_itrs.size(); ++i) {
                color_id_itrs[i].next_geq(candidate);
                color_t val = color_id_itrs.at(i).value();
                if (val != candidate) {
                    candidate = val;
                    i = 0;
                    break;
                }
            }
            if (i == color_id_itrs.size()) {
                colors.push_back(candidate);
                color_id_itrs[0].next();
                candidate = color_id_itrs.at(0).value();
                i = 1;
            }
        }
    }
    return colors;
}


#undef CLASS_HEADER
#undef METHOD_HEADER

} // namespace minimizer
} // namespace kaminari