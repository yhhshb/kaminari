#include "../../include/index.hpp"
#include "../../include/hybrid.hpp"


namespace kaminari {
namespace minimizer {

#define CLASS_HEADER template <class ColorClasses, class ColorMapper>
#define METHOD_HEADER index<ColorClasses, ColorMapper>

template class index<kaminari::color_classes::hybrid, pthash::compact_vector>;

auto bit_to_nuc(uint64_t bits, uint32_t len) {
    std::string nuc;
    for (uint32_t i = 0; i < len; ++i) {
        uint8_t b = (bits >> (2 * (len - i - 1))) & 3;
        switch (b) {
            case 0: nuc += 'A'; break;
            case 1: nuc += 'C'; break;
            case 2: nuc += 'G'; break;
            case 3: nuc += 'T'; break;
        }
    }
    return nuc;
}


CLASS_HEADER
std::vector<scored_id> 
METHOD_HEADER::ranking_query_union_threshold(char const * const q, const std::size_t l, float threshold_ratio, std::size_t verbosity_level) const noexcept
{
    if (verbosity_level > 1) std::cerr << "step 1: collect color class ids\n";
    if (verbosity_level > 2) std::cerr << "s : " << q << " l : " << l << "\n"; 
    std::vector<std::pair<std::size_t, uint32_t>> ccids_counts;
    uint64_t contig_kmer_count; 
    { // collect color class ids
        std::size_t contig_mmer_count;
        std::vector<::minimizer::record_t> mms_buffer;
        contig_kmer_count = ::minimizer::from_string<hash64>(q, l, k, m, seed, canonical, contig_mmer_count, mms_buffer);
        for (const auto& record : mms_buffer) { 
            //std::cerr << "record : " << record.itself << " -> " << bit_to_nuc(record.itself, m) << " -pthash> " <<  hf(record.itself) <<  " -ccid> " <<  m_map[hf(record.itself)] << " -size> " << record.size << "\n";
            ccids_counts.push_back(std::make_pair(m_map[hf(record.itself)], record.size));
        }
        if (verbosity_level > 3) std::cerr << "query contains " << contig_kmer_count << " k-mers and " << contig_mmer_count << " m-mers\n";
    }
    if (verbosity_level > 3) std::cerr << ccids_counts << "\n";


    /* std::vector<uint32_t> pbs;
    uint64_t tot = 0;
    for (int i = 0; i<ccids_counts.size(); i++){
        std::cerr << ccids_counts[i].first << " " << ccids_counts[i].second << "\n";
        typename ColorClasses::row_accessor row = m_ccs.colors_at(ccids_counts[i].first);
        while(row.value() != m_filenames.size()){
            //std::cerr << row.value() << "\n";
            if (row.value() == 914){
                pbs.push_back(i);
                tot += ccids_counts[i].second;
            } 
            row.next();
        }
        row.reset();        
    }
    std::cerr << "pbs : " << pbs << "\n";
    std::cerr << "tot : " << tot << "\n"; */


    if (verbosity_level > 2) std::cerr << "step 2: ids to colors\n";
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

        if (verbosity_level > 3) std::cerr << ccids_counts << " (post sort)\n";

        
    

        for (auto itr = ccids_counts.begin(); itr != last; ++itr) { 
            color_itrs.push_back(std::make_pair(m_ccs.colors_at((*itr).first), (*itr).second));
            if (color_itrs.back().first.type() != ColorClasses::row_accessor::complementary_delta_gaps) {
                all_very_dense = false;
            }
        }
    }

    
    
    if (verbosity_level > 2) {
        std::cerr << "step 3: computing intersections\n";
        if (all_very_dense) std::cerr << "\tcompute dense intersection\n";
        else std::cerr << "\tcompute mixed intersection (for a mix of dense and sparse vectors)\n"; 
    }
    if (color_itrs.empty()) return {};
    //if (all_very_dense) return ranking_dense_intersection(std::move(color_itrs), contig_kmer_count*threshold_ratio, contig_kmer_count, verbosity_level); // intersect of dense rows
    return ranking_mixed_intersection(std::move(color_itrs), contig_kmer_count*threshold_ratio, verbosity_level); // intersect dense and sparse rows
}


CLASS_HEADER
std::vector<typename METHOD_HEADER::color_t> 
METHOD_HEADER::ranking_dense_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold, uint64_t nb_kmers, std::size_t verbosity_level) const noexcept
{
    for (auto i : color_id_itrs)
    {
        std::cerr << i.first.size() << "\n";
    }

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
std::vector<scored_id> 
METHOD_HEADER::ranking_mixed_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold, std::size_t verbosity_level) const noexcept
{
    if (threshold == 0) threshold = 1; //super low values of opts.threshold_ratio

    std::sort(
        color_id_itrs.begin(),
        color_id_itrs.end(),
        [](auto const& x, auto const& y) { return x.first.size() < y.first.size(); }
    );


    /* int ind = 0;
    std::vector<uint32_t> pbs;
    uint64_t tot = 0;
    for (auto i : color_id_itrs)
    {
        std::cerr << i.first.size() << " size, " << i.second << " nb kmers \n";
        while(i.first.value() != m_filenames.size()){
            if (i.first.value() == 0){
                pbs.push_back(ind);
                tot += i.second;
            }
            //if (i.first.value() == 0) std::cerr << "0 found\n";
            //std::cerr << i.first.value() << "\n";
            i.first.next();
        }
        i.first.reset();     
        ind ++;   
    }
    std::cerr << "pbs : " << pbs << "\n";
    std::cerr << "tot : " << tot << "\n"; */


    bit::vector<uint64_t> tested(m_filenames.size());
    std::vector<scored_id> colors;
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
                    candidate = color_id_itrs[i].first.value();
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
                    color_t val = color_id_itrs[idx].first.value();
                    if (val == candidate) {
                        score += color_id_itrs[idx].second;
                    } 
                }
                
                if (score >= threshold) {
                    colors.push_back({candidate, score});
                }

                color_id_itrs[i].first.next();
                candidate = color_id_itrs[i].first.value();
            }
        }
    }
    
    return colors;
}


#undef CLASS_HEADER
#undef METHOD_HEADER

} // namespace minimizer
} // namespace kaminari