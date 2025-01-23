#include "../../include/index.hpp"
#include "../../include/hybrid.hpp"



namespace kaminari {
namespace minimizer {

#define CLASS_HEADER template <class ColorClasses, class ColorMapper>
#define METHOD_HEADER index<ColorClasses, ColorMapper>

template class index<kaminari::color_classes::hybrid, pthash::compact_vector>;

CLASS_HEADER
std::vector<typename METHOD_HEADER::color_t> 
METHOD_HEADER::query_union_threshold(char const * const q, const std::size_t l, options_t& opts) const noexcept
{
    if (opts.verbose > 1) std::cerr << "step 1: collect color class ids\n";
    if (opts.verbose > 2) std::cerr << "s : " << q << " l : " << l << "\n"; 
    std::vector<std::pair<std::size_t, uint32_t>> ccids_counts;
    uint64_t contig_kmer_count; 
    { // collect color class ids
        std::size_t contig_mmer_count;
        std::vector<::minimizer::record_t> mms_buffer;
        contig_kmer_count = ::minimizer::from_string<hash64>(q, l, k, m, seed, canonical, contig_mmer_count, mms_buffer);
        for (const auto& record : mms_buffer) { 
            //std::cerr << "record : " << record.itself << "\n";
            uint32_t cid_with_parity = m_map[hf(record.itself)];
            if ((record.itself & 1) == (cid_with_parity & 1)){
                //checkin not alien kmer
                ccids_counts.push_back(std::make_pair(cid_with_parity >> 1, record.size)); //masking out parity
            }
        }
        if (opts.verbose > 3) std::cerr << "query contains " << contig_kmer_count << " k-mers and " << contig_mmer_count << " m-mers\n";
    }
    if (opts.verbose > 3) std::cerr << ccids_counts << "\n";
    if (opts.verbose > 2) std::cerr << "step 2: ids to colors\n";
    std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>> color_itrs;
    bool all_very_dense = true;


    {
        //sort and unique+erase to remove duplicates, side effect : accumulate counts (n°kmers) for every duplicate
        std::sort(
            ccids_counts.begin(),
            ccids_counts.end(),
            [](auto const& x, auto const& y) { return x.first < y.first; }
        );
        auto last = ccids_counts.erase( kaminari::utils::unique_accumulate( 
                                            ccids_counts.begin(), 
                                            ccids_counts.end() ), 
                                        ccids_counts.end() );

        if (opts.verbose > 3) std::cerr << ccids_counts << " (post sort)\n";

        for (auto itr = ccids_counts.begin(); itr != last; ++itr) { 
            color_itrs.push_back(std::make_pair(m_ccs.colors_at((*itr).first), (*itr).second));
            if (color_itrs.back().first.type() != ColorClasses::row_accessor::complementary_delta_gaps) {
                all_very_dense = false;
            }
        }
    }
    
    if (opts.verbose > 2) { 
        std::cerr << "step 3: computing intersections\n";
        if (all_very_dense) std::cerr << "\tcompute dense intersection\n";
        else std::cerr << "\tcompute mixed intersection (for a mix of dense and sparse vectors)\n"; 
    }
    if (color_itrs.empty()) return {};
    if (all_very_dense) return union_dense_intersection(std::move(color_itrs), contig_kmer_count*opts.threshold_ratio, contig_kmer_count); // intersect of dense rows
    
    return union_mixed_intersection(std::move(color_itrs), contig_kmer_count*opts.threshold_ratio); // intersect dense and sparse rows
}


CLASS_HEADER
std::vector<typename METHOD_HEADER::color_t> 
METHOD_HEADER::union_dense_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold, uint64_t nb_kmers) const noexcept
        /* check every complementary (= absent docids from color) and for each color
         where it is absent, remove color_count to its count (because the candidate is absent to this number of kmers)
         at the end, if the candidate has been absent to too many kmers (threshold), it is not chosen
        */
{
    if (threshold == 0) threshold = 1; //super low values of opts.threshold_ratio
    
    std::vector<color_t> colors;
    std::size_t vec_size = color_id_itrs.size();
    std::size_t filenames_size = m_filenames.size();

    std::vector<uint32_t> counts(filenames_size, 0);
    uint64_t global_count = 0;

    for (uint64_t i = 0; i != vec_size; ++i) {
        global_count += color_id_itrs[i].second;
        while (color_id_itrs[i].first.comp_value() < filenames_size) { //"<" and not "!=" because in case of color having every docid, comp_val stays -1, in an unsigned int so 2^32-1
            counts[color_id_itrs[i].first.comp_value()] += color_id_itrs[i].second;
            color_id_itrs[i].first.comp_next();
        }
    }

    for (uint32_t i = 0; i != filenames_size; ++i) {
        if (global_count - counts[i] >= threshold) colors.push_back(i);
    }

    return colors;
}


CLASS_HEADER
std::vector<typename METHOD_HEADER::color_t>  
METHOD_HEADER::union_mixed_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold) const noexcept
{
    if (threshold == 0) threshold = 1; //super low values of opts.threshold_ratio
    
    std::vector<color_t> colors;
    std::size_t vec_size = color_id_itrs.size();
    std::size_t filenames_size = m_filenames.size();

    std::vector<uint32_t> counts(filenames_size, 0);

    for (uint64_t i = 0; i != vec_size; ++i) {
        while (color_id_itrs[i].first.value() != filenames_size) {
            counts[color_id_itrs[i].first.value()] += color_id_itrs[i].second;
            color_id_itrs[i].first.next();
        }
    }

    for (uint32_t i = 0; i != filenames_size; ++i) {
        if (counts[i] >= threshold) colors.push_back(i);
    }

    return colors;
} 


#undef CLASS_HEADER
#undef METHOD_HEADER

} // namespace minimizer
} // namespace kaminari
