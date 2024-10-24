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
            ccids_counts.push_back(std::make_pair(m_map[hf(record.itself)], record.size));
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
    //if (all_very_dense) return union_dense_intersection(std::move(color_itrs), contig_kmer_count*opts.threshold_ratio, contig_kmer_count, opts.verbose); // intersect of dense rows

    //TODO : opti dense case aswell, for now only mixed has been optimized so use it for experiments
    return union_mixed_intersection(std::move(color_itrs), contig_kmer_count*opts.threshold_ratio, opts.verbose); // intersect dense and sparse rows
}


CLASS_HEADER
std::vector<typename METHOD_HEADER::color_t> 
METHOD_HEADER::union_dense_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold, uint64_t nb_kmers, std::size_t verbosity_level) const noexcept
{
    if (threshold == 0) threshold = 1;
    std::size_t filenames_size = m_filenames.size();
    std::vector<color_t> tmp;

    { // step 1: take the union of complementary sets
        for (auto& itr : color_id_itrs) itr.first.reinit_for_complemented_set_iteration();
        color_t candidate = (
            *std::min_element(
                color_id_itrs.begin(), 
                color_id_itrs.end(),
                [](auto const& x, auto const& y) {
                    return x.first.comp_value() < y.first.comp_value();
                }
            )
        ).first.comp_value();

        /* check every complementary (= absent docids from color) and for each color
         where it is absent, remove color_count to its count (because the candidate is absent to this number of kmers)
         at the end, if the candidate has been absent to too many kmers (threshold), it is put in tmp (= eliminated)
         everything that has not been inserted in tmp, has else a sufficient count, or hasn't been a complementary == present everywhere) 
        */
        uint32_t candidate_count;
        tmp.reserve(filenames_size);
        while (candidate < filenames_size) {
            candidate_count = nb_kmers;
            color_t next_candidate = filenames_size;
            for (uint64_t i = 0; i != color_id_itrs.size(); ++i) {
                if (color_id_itrs.at(i).first.comp_value() == candidate){
                    candidate_count -= color_id_itrs.at(i).second;
                    color_id_itrs[i].first.comp_next();
                } 
                // compute next minimum 
                if (color_id_itrs.at(i).first.comp_value() < next_candidate) {
                    next_candidate = color_id_itrs.at(i).first.comp_value();
                }
            }
            if (candidate_count < threshold) tmp.push_back(candidate);
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
METHOD_HEADER::union_mixed_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold, std::size_t verbosity_level) const noexcept
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

    for (uint64_t i = 0; i != filenames_size; ++i) {
        if (counts[i] >= threshold) colors.push_back(i);
    }

    return colors;
} 


/* FULGOR's way
CLASS_HEADER
std::vector<typename METHOD_HEADER::color_t>  
METHOD_HEADER::union_mixed_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold, std::size_t verbosity_level) const noexcept
{
    if (threshold == 0) threshold = 1; //super low values of opts.threshold_ratio

    std::sort(
        color_id_itrs.begin(),
        color_id_itrs.end(),
        [](auto const& x, auto const& y) { return x.first.size() < y.first.size(); }
    );
    
    std::vector<color_t> colors;
    std::size_t vec_size = color_id_itrs.size();
    std::size_t filenames_size = m_filenames.size();
    std::size_t score;

    color_t candidate =
        (*std::min_element(color_id_itrs.begin(), color_id_itrs.end(), [](auto const& x, auto const& y) {
            return x.first.value() < y.first.value();
        })).first.value();

    uint64_t nb_candidates = 0;

    while (candidate != filenames_size) {
        nb_candidates++;
        uint32_t next_candidate = filenames_size;
        uint32_t score = 0;
        for (uint64_t i = 0; i != vec_size; ++i) {
            if (color_id_itrs[i].first.value() == candidate) {
                score += color_id_itrs[i].second;
                color_id_itrs[i].first.next();
            }
            // compute next minimum
            if (color_id_itrs[i].first.value() < next_candidate) {
                next_candidate = color_id_itrs[i].first.value();
            }
        }
        if (score >= threshold) colors.push_back(candidate);
        assert(next_candidate > candidate);
        candidate = next_candidate;
    }

    return colors;
} */



/* V1 : heuristical way 
CLASS_HEADER
std::vector<typename METHOD_HEADER::color_t>  
METHOD_HEADER::union_mixed_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold, std::size_t verbosity_level) const noexcept
{
    if (threshold == 0) threshold = 1; //super low values of opts.threshold_ratio

    std::sort(
        color_id_itrs.begin(),
        color_id_itrs.end(),
        [](auto const& x, auto const& y) { return x.first.size() < y.first.size(); }
    );

    //for (auto i : color_id_itrs)
    {
        std::cerr << "next i \n";
        std::cerr << i.first.size() << " size, " << i.second << " nb kmers \n";
        while(i.first.value() != m_filenames.size()){
            //std::cerr << i.first.value() << "\n";
            i.first.next();
        }
        i.first.reset();        
    }
    

    bit::vector<uint64_t> tested(m_filenames.size());
    std::vector<color_t> colors;
    std::size_t vec_size = color_id_itrs.size();
    std::size_t filenames_size = m_filenames.size();
    std::size_t idx;
    std::size_t score;
    bit::ef::array remaining_counts;

    {
        //used to know if worth looking in next colors
        std::vector<uint64_t> ef_build(vec_size+1, 0);
        uint64_t sum = 0;
        for (int i = 1; i <= vec_size; i++) {
            sum += color_id_itrs[vec_size-i].second;
            ef_build[i] += sum;
        }
        ef_build[vec_size] += 1; //acts like a boundary
        remaining_counts = bit::ef::array(ef_build.begin(), vec_size+1, ef_build.back());
    }

    {
        color_t candidate;
        //stop loop when not enough kmers represented in following row_acc
        //to reach threshold 
        for (size_t i = 0; i <= vec_size-1 - remaining_counts.lt_find(threshold); i++){
            if (i>0) color_id_itrs[i].first.reset();
            candidate = color_id_itrs[i].first.value();

            //stop loop when regularly reaching last element
            while (candidate != filenames_size){

                // check if candidate has not been done in other previous row_acc
                //+ need to check if not filenames_size because of tested out_of_range 
                while (candidate != filenames_size && tested.at(candidate)) {
                    color_id_itrs[i].first.next();
                    candidate = color_id_itrs.at(i).first.value();
                }
                // leaving loop when heuristically reaching last element
                if (candidate == filenames_size) {
                    break;
                    //checked every value in this row accessor before next i
                }
                     
                tested.set(candidate); //we chose a never seen candidate, dont check him again

                idx = i+1;
                score = color_id_itrs[i].second;
                
                //for (; (idx < vec_size) && (score < threshold); ++idx) { (removed heuristics to stop when score reached threshol for expe FP diff to truth)
                for (; (idx < vec_size); ++idx) {
                    if (i>0) color_id_itrs[idx].first.reset(); //candidate in 2nd row_acc might be lower than candidates in 1st row_acc
                    color_id_itrs[idx].first.next_geq(candidate);
                    color_t val = color_id_itrs.at(idx).first.value();
                    if (val == candidate) {
                        score += color_id_itrs.at(idx).second;
                    } 
                }
                
                if (score >= threshold) {
                    colors.push_back(candidate);
                }

                color_id_itrs.at(i).first.next();
                candidate = color_id_itrs[i].first.value();
            }
        }
    }
    
    return colors;
} */


#undef CLASS_HEADER
#undef METHOD_HEADER

} // namespace minimizer
} // namespace kaminari