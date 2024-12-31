#include "../../include/index.hpp"
#include "../../include/hybrid.hpp"

namespace kaminari {
namespace minimizer {

#define CLASS_HEADER template <class ColorClasses, class ColorMapper>
#define METHOD_HEADER index<ColorClasses, ColorMapper>

template class index<kaminari::color_classes::hybrid, pthash::compact_vector>;


// Function to read a file
CLASS_HEADER
void 
METHOD_HEADER::read_file_task(const std::string& file, uint32_t doc_id, Depth1Result& result, const build::options_t& build_parameters, std::mutex& debug_cerr_mutex) {
    std::size_t total_kmers = 0;
    std::size_t total_mmers = 0;
    std::size_t total_minimizers = 0;
    gzFile fp = nullptr;
    kseq_t* seq = nullptr;
    std::vector<::minimizer::record_t> mms_buffer;
    if ((fp = gzopen(file.c_str(), "r")) == NULL)
        throw std::runtime_error("Unable to open input file " + file);
    seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        std::size_t contig_mmer_count;
        auto contig_kmer_count = ::minimizer::from_string<hash64>(
            seq->seq.s, 
            seq->seq.l, 
            build_parameters.k, 
            build_parameters.m, 
            build_parameters.seed, 
            build_parameters.canonical,
            contig_mmer_count,
            mms_buffer
        );
        total_kmers += contig_kmer_count;
        total_mmers += contig_mmer_count;
        total_minimizers += mms_buffer.size();
        if (build_parameters.verbose > 4) {
            std::lock_guard<std::mutex> lock(debug_cerr_mutex);
            std::cerr << "\tRead contig:\n" 
                        << "\t\t" << contig_kmer_count << " k-mers\n"
                        << "\t\t" << contig_mmer_count << " m-mers\n"
                        << "\t\t" << mms_buffer.size() << " minimizers\n";
        }
        // duplicates removed during merge
        std::vector<minimizer_t> minimizers;
        for (auto r : mms_buffer) {
            result.minmers.push_back(r.itself);
        }
        mms_buffer.clear();
    }
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    fp = nullptr;
    seq = nullptr;
    result.docid = doc_id;
}



// Function to merge two or 3 depth_1 results
CLASS_HEADER
void
METHOD_HEADER::merge_results_task(const Depth1Result& left, const Depth1Result& right, depthn_result& result, std::mutex& debug_cerr_mutex) {
    // Determine insertion order based on docid value
    const Depth1Result* small = &left;
    const Depth1Result* big = &right;
    if (left.docid > right.docid) {
        small = &right;
        big = &left;
    }

    auto small_itr = small->minmers.cbegin();  
    auto big_itr = big->minmers.cbegin(); 
    uint64_t last = 1ULL<<63;
    //used for deduplicating, first ensuring last isnt equal to first value considered

    //read both vectors
    while (small_itr != small->minmers.cend() && big_itr != big->minmers.cend()) {
        if (*small_itr < *big_itr) {
            if (*small_itr != last) {
                result.push_back(std::make_pair(*small_itr, std::vector<uint32_t> {small->docid}));
                last = *small_itr;
            }
            ++small_itr;
        } else if (*small_itr > *big_itr) {
            if (*big_itr != last) {
                result.push_back(std::make_pair(*big_itr, std::vector<uint32_t> {big->docid}));
                last = *big_itr;
            }
            ++big_itr;
        } else {
            if (*small_itr != last) {
                result.push_back(std::make_pair(*small_itr, std::vector<uint32_t> {small->docid, big->docid}));
                last = *small_itr;
            }
            ++small_itr;
            ++big_itr;
        }
    } 

    //read remaining non-empty vector
    while (small_itr != small->minmers.cend()) {
        if (*small_itr != last) {
            result.push_back(std::make_pair(*small_itr, std::vector<uint32_t> {small->docid}));
            last = *small_itr;
        }
        ++small_itr;
    }

    while (big_itr != big->minmers.cend()) {
        if (*big_itr != last) {
            result.push_back(std::make_pair(*big_itr, std::vector<uint32_t> {big->docid}));
            last = *big_itr;
        }
        ++big_itr;
    }
}
CLASS_HEADER
void
METHOD_HEADER::merge_results_task(const Depth1Result& one, const Depth1Result& two, Depth1Result& three, depthn_result& result, std::mutex& debug_cerr_mutex) {
    const Depth1Result* first = &one;
    const Depth1Result* second = &two;
    const Depth1Result* third = &three;

    if (first->docid > second->docid) std::swap(first, second);
    if (first->docid > third->docid) std::swap(first, third);
    if (second->docid > third->docid) std::swap(second, third);

    auto first_itr = first->minmers.cbegin();
    auto second_itr = second->minmers.cbegin();
    auto third_itr = third->minmers.cbegin();

    uint64_t last = 1ULL<<63; // Used for deduplication; initialize to an impossible value

    // Read all three vectors
    while (first_itr != first->minmers.cend() || second_itr != second->minmers.cend() || third_itr != third->minmers.cend()) {
        uint64_t min_value = 1ULL<<63;

        //min between three
        if (first_itr != first->minmers.cend()) min_value = std::min(min_value, *first_itr);
        if (second_itr != second->minmers.cend()) min_value = std::min(min_value, *second_itr);
        if (third_itr != third->minmers.cend()) min_value = std::min(min_value, *third_itr);

        std::vector<uint32_t> doc_ids;

        //collect docids associated with min value
        if (first_itr != first->minmers.cend() && *first_itr == min_value) {
            doc_ids.push_back(first->docid);
            ++first_itr;
        }
        if (second_itr != second->minmers.cend() && *second_itr == min_value) {
            doc_ids.push_back(second->docid);
            ++second_itr;
        }
        if (third_itr != third->minmers.cend() && *third_itr == min_value) {
            doc_ids.push_back(third->docid);
            ++third_itr;
        }
        // deduplicate if needed
        if (min_value != last) {
            result.push_back(std::make_pair(min_value, std::move(doc_ids)));
            last = min_value;
        }
    }
}

//special case of merging depth1 result with depthn result (needed when 1 file in last batch)
CLASS_HEADER
void
METHOD_HEADER::merge_results_task(const depthn_result& left, const Depth1Result& right, depthn_result& result, std::mutex& debug_cerr_mutex) {
    auto left_itr = left.cbegin();     // depthn_result is a deque of <vector<uint32_t>, uint64_t>
    auto right_itr = right.minmers.cbegin();  // Depth1Result's `first` is a vector of uint64_t

    uint64_t right_docid = right.docid;   // Depth1Result's `second` is the document ID
    uint64_t last = 1ULL<<63; // Deduplication marker

    // Merge logic
    while (right_itr != right.minmers.cend() || left_itr != left.cend()) {
        uint64_t min_value = 1ULL<<63;
        if (left_itr != left.cend()) min_value = (*left_itr).first;
        if (right_itr != right.minmers.cend() && *right_itr < min_value) min_value = *right_itr;
        
        std::vector<uint32_t> merged_docids;

        // Add from left if it matches the smallest value
        if (left_itr != left.cend() && (*left_itr).first == min_value) {
            merged_docids.insert(merged_docids.end(), (*left_itr).second.begin(), (*left_itr).second.end());
            ++left_itr;
        }

        // Add from right if it matches the smallest value
        if (right_itr != right.minmers.cend() && *right_itr == min_value) {
            merged_docids.push_back(right_docid);
            ++right_itr;
        }

        // Deduplicate and insert into the result
        if (min_value != last) {
            result.push_back(std::make_pair(min_value, std::move(merged_docids)));
            last = min_value;
        }
    }
}


// Functions to merge two or 3 Depth n results
CLASS_HEADER
void
METHOD_HEADER::merge_results_task(const depthn_result& left, const depthn_result& right, depthn_result& result, std::mutex& debug_cerr_mutex) {
    // Determine insertion order based on docid value
    const depthn_result* small = &left;
    const depthn_result* big = &right;
    //check 1st value of each first color
    if ((*left.cbegin()).second[0] > (*right.cbegin()).second[0]) {
        small = &right;
        big = &left;
    }

    // Merge the two results
    auto small_itr = small->cbegin();  
    auto big_itr = big->cbegin(); 

    while (small_itr != small->cend() && big_itr != big->cend()) {
        if ((*small_itr).first < (*big_itr).first) {
            result.push_back(*small_itr);
            ++small_itr;
        } else if ((*small_itr).first > (*big_itr).first) {
            result.push_back(*big_itr);
            ++big_itr;
        } else {
            std::vector<uint32_t> merged_docids = (*small_itr).second;
            merged_docids.insert(merged_docids.end(), (*big_itr).second.begin(), (*big_itr).second.end()); // merging by concatenating
            result.push_back(std::make_pair((*small_itr).first, std::move(merged_docids))); // Add to result
            ++small_itr;
            ++big_itr;
        }
    }

    while (small_itr != small->cend()) {
        result.push_back(*small_itr);
        ++small_itr;
    }
    while (big_itr != big->cend()) {
        result.push_back(*big_itr);
        ++big_itr;
    }
}
CLASS_HEADER
void
METHOD_HEADER::merge_results_task(const depthn_result& one, const depthn_result& two, const depthn_result& three, depthn_result& result, std::mutex& debug_cerr_mutex) {
    // Determine insertion order based on docid value
    const depthn_result* first = &one;
    const depthn_result* second = &two;
    const depthn_result* third = &three;

    //check order based on first value of first color
    if ((*second->cbegin()).second[0] < (*first->cbegin()).second[0]) {
        std::swap(first, second);
    }
    if ((*third->cbegin()).second[0] < (*first->cbegin()).second[0]) {
        std::swap(first, third);
    }
    if ((*third->cbegin()).second[0] < (*second->cbegin()).second[0]) {
        std::swap(second, third);
    }

    // Iterators for the three inputs
    auto itr1 = first->cbegin();
    auto itr2 = second->cbegin();
    auto itr3 = third->cbegin();

    // Merge the three results
    while (itr1 != first->cend() || itr2 != second->cend() || itr3 != third->cend()) {
        uint64_t min_value = 1ULL<<63;
        if (itr1 != first->cend()) min_value = (*itr1).first;
        if (itr2 != second->cend() && (*itr2).first < min_value) min_value = (*itr2).first;
        if (itr3 != third->cend() && (*itr3).first < min_value) min_value = (*itr3).first;

        std::vector<uint32_t> merged_docids;

        if (itr1 != first->cend() && (*itr1).first == min_value) {
            merged_docids.insert(merged_docids.end(), (*itr1).second.begin(), (*itr1).second.end());
            ++itr1;
        }
        if (itr2 != second->cend() && (*itr2).first == min_value) {
            merged_docids.insert(merged_docids.end(), (*itr2).second.begin(), (*itr2).second.end());
            ++itr2;
        }
        if (itr3 != third->cend() && (*itr3).first == min_value) {
            merged_docids.insert(merged_docids.end(), (*itr3).second.begin(), (*itr3).second.end());
            ++itr3;
        }

        // Add the merged result to the output
        result.push_back(std::make_pair(min_value, std::move(merged_docids)));
    }
}

//function to merge 2 depthn results into final result
CLASS_HEADER
void
METHOD_HEADER::merge_results_task(const depthn_result& left, const depthn_result& right, colors_to_minmer& final_result, std::mutex& debug_cerr_mutex) {
    // Determine insertion order based on docid value
    const depthn_result* small = &left;
    const depthn_result* big = &right;
    //check 1st value of each first color
    if ((*left.cbegin()).second[0] > (*right.cbegin()).second[0]) {
        small = &right;
        big = &left;
    }

    // Merge the two results
    auto small_itr = small->cbegin();  
    auto big_itr = big->cbegin(); 

    while (small_itr != small->cend() && big_itr != big->cend()) {
        if ((*small_itr).first < (*big_itr).first) {
            final_result.push_back(std::make_pair((*small_itr).second, (*small_itr).first));
            ++small_itr;
        } else if ((*small_itr).first > (*big_itr).first) {
            final_result.push_back(std::make_pair((*big_itr).second, (*big_itr).first));
            ++big_itr;
        } else {
            std::vector<uint32_t> merged_docids = (*small_itr).second;
            merged_docids.insert(merged_docids.end(), (*big_itr).second.begin(), (*big_itr).second.end()); // merging by concatenating
            final_result.push_back(std::make_pair(std::move(merged_docids), (*small_itr).first)); // Add to final_result
            ++small_itr;
            ++big_itr;
        }
    }
    while (small_itr != small->cend()) {
        final_result.push_back(std::make_pair((*small_itr).second, (*small_itr).first));
        ++small_itr;
    }
    while (big_itr != big->cend()) {
        final_result.push_back(std::make_pair((*big_itr).second, (*big_itr).first));
        ++big_itr;
    }
}


#undef CLASS_HEADER
#undef METHOD_HEADER

} // namespace minimizer
} // namespace kaminari