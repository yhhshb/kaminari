#ifndef RANKING_THRESHOLD_UNION_QUERY_HPP
#define RANKING_THRESHOLD_UNION_QUERY_HPP

#include "../utils.hpp"

namespace kaminari {
namespace minimizer {

std::vector<scored_id> 
index::ranking_query_union_threshold(colorsets_accessor& parser, char const * const q, const std::size_t l, options_t& opts) const noexcept
{
    if (opts.verbose > 1) std::cerr << "step 1: collect color class ids\n";
    if (opts.verbose > 2) std::cerr << "s : " << q << " l : " << l << "\n"; 
    std::vector<std::pair<std::uint32_t, uint32_t>> ccids_counts;
    uint64_t contig_kmer_count; 

    { // collect color class ids
        std::size_t contig_mmer_count;
        std::vector<::minimizer::record_t> mms_buffer;
        contig_kmer_count = ::minimizer::from_string<hash64>(q, l, k, m, seed, canonical, contig_mmer_count, mms_buffer);
        for (const auto& record : mms_buffer) { 
            uint32_t cid_with_parity = m_colormapper[hf(record.hash)];
            //std::cerr << record.hash << " -> " << record.size << " -> " << hf(record.hash) << " -> " << cid_with_parity << "\n";
            if ((record.hash & ((1UL << b)-1)) == (cid_with_parity & ((1UL << b)-1))){
                //checkin not alien kmer
                ccids_counts.push_back(std::make_pair(cid_with_parity >> b, record.size)); //masking out parity
            }
        }


        if (opts.verbose > 3) std::cerr << "query contains " << contig_kmer_count << " k-mers and " << contig_mmer_count << " m-mers\n";
    }
    if (opts.verbose > 3) std::cerr << ccids_counts << "\n";

    if (opts.verbose > 2) std::cerr << "step 2: ids to colors\n";


    {
        //sort and unique+erase to remove duplicates, side effect : accumulate counts (n°kmers) for every duplicate
        std::sort(
            ccids_counts.begin(),
            ccids_counts.end(),
            [](auto const& x, auto const& y) { return x.first < y.first; }
        );
        ccids_counts.erase( 
            kaminari::utils::unique_accumulate( 
                ccids_counts.begin(), 
                ccids_counts.end() ), 
            ccids_counts.end() );

        if (opts.verbose > 3) std::cerr << ccids_counts << " (post sort)\n";
    }

    if (opts.verbose > 2) {
        std::cerr << "step 3: computing intersections\n";
    }
    if (ccids_counts.empty()) return {};
    return ranking_mixed_intersection(parser, std::move(ccids_counts), contig_kmer_count*opts.threshold_ratio); // intersect dense and sparse rows
}


std::vector<scored_id> 
index::ranking_mixed_intersection(colorsets_accessor& parser,std::vector<std::pair<std::uint32_t, uint32_t>>&& ccids_counts, uint64_t threshold) const noexcept
{
    if (threshold == 0) threshold = 1; //super low values of opts.threshold_ratio
    
    std::vector<scored_id> colors;
    std::size_t vec_size = ccids_counts.size();

    std::vector<uint32_t> counts(nb_docs, 0);
    
    for (size_t i = 0; i != vec_size; ++i) {
        parser.move_to(ccids_counts[i].first);
        while (parser.value() != nb_docs) {
            counts[parser.value()] += ccids_counts[i].second;
            parser.next();
        }
    }

    for (uint32_t i = 0; i != nb_docs; ++i) {
        if (counts[i] >= threshold) colors.push_back(scored_id{i, counts[i]});
    }

    return colors;
    
}

} // namespace minimizer
} // namespace kaminari

#endif // RANKING_THRESHOLD_UNION_QUERY_HPP
