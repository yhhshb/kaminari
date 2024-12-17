#ifndef KAMINARI_INDEX_HPP
#define KAMINARI_INDEX_HPP

#include <iostream>
#include <mutex>
#include <shared_mutex>
#include <condition_variable>
#include <deque>
#include <stack>
#include "constants.hpp"
#include "utils.hpp"
#include "build_options.hpp"
#include "query_options.hpp"
#include "minimizer.hpp"
#include "../bundled/pthash/include/pthash.hpp"
#include "../bundled/biolib/bundled/prettyprint.hpp"
#include "../bundled/biolib/include/bit_vector.hpp"
#include "../bundled/biolib/include/elias_fano.hpp"
#include "../bundled/unordered_dense/include/ankerl/unordered_dense.h"

#include <chrono>

namespace kaminari {
namespace minimizer {


#define CLASS_HEADER template <class ColorClasses, class ColorMapper>
#define METHOD_HEADER index<ColorClasses, ColorMapper>


template <typename T>
struct scored {
    T item;
    uint32_t score;
};

typedef scored<uint32_t> scored_id;

typedef emem::external_memory_vector<uint64_t> emem_minmers;
struct Depth1Result {
    emem_minmers minmers;
    uint32_t docid;
};

typedef emem::external_memory_vector<std::pair<uint64_t, std::vector<uint32_t>>, false> depthn_result;

typedef emem::external_memory_vector<std::pair<std::vector<uint32_t>, uint64_t>> colors_to_minmer;

struct Function {
    std::function<void()> f; // Actual function
    std::string desc;
};

// Custom hash function for vector<uint32_t>
struct VectorHash {
    std::size_t operator()(std::vector<uint32_t> const& vec) const {
        std::size_t seed = vec.size();
        for(auto x : vec) {
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = (x >> 16) ^ x;
            seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

CLASS_HEADER
class index
{
    public:
        typedef uint64_t minimizer_t;
        typedef typename ColorClasses::color_t color_t;
        typedef typename kaminari::query::options_t options_t;

        index();
        index(const build::options_t& build_parameters);

        //following methods are explicitly instantiated in src/psa/
        std::vector<color_t> query_full_intersection(char const * const q, const std::size_t l, options_t& opts) const noexcept;

        std::vector<color_t> query_union_threshold(char const * const q, const std::size_t l, options_t& opts) const noexcept;

        std::vector<scored_id> ranking_query_union_threshold(char const * const q, const std::size_t l, options_t& opts) const noexcept;
        
        void memory_breakdown(std::ostream& out) const noexcept;
        
        template <class Visitor>
        void visit(Visitor& visitor);

    private:
        typedef pthash::build_configuration pthash_opt_t;
        typedef pthash::dense_partitioned_phf<pthash::murmurhash2_64, pthash::opt_bucketer, pthash::mono_EF, true, pthash::pthash_search_type::add_displacement> pthash_minimizers_mphf_t;

        template <class Iterator>
        class pthash_input_iterator {
            public:
                pthash_input_iterator(Iterator mm_itr) : m_iterator(mm_itr) {}
                void operator++() {++m_iterator;}
                typename Iterator::value_type operator*() const {return (*m_iterator);}
            private:
                Iterator m_iterator;
        };

        pthash_opt_t get_pthash_options(const build::options_t& build_parameters);
        void build(const build::options_t& build_parameters);
        void build2(const build::options_t& build_parameters);
        void worker_thread(
            std::stack<Function>& task_stack,
            std::mutex& task_stack_mutex,
            std::condition_variable& cv,
            std::atomic<uint32_t>& running_task,
            std::atomic<bool>& all_done,
            std::mutex& debug_cerr_mutex);

        //merge 2 or 3 (case odd number) results post parsing files
        void merge_results_task(const Depth1Result& left, const Depth1Result& right, depthn_result& result, std::mutex& debug_cerr_mutex);
        void merge_results_task(const Depth1Result& one, const Depth1Result& two, Depth1Result& three, depthn_result& result, std::mutex& debug_cerr_mutex);
        //merge 2 or 3 (case odd number) results 
        void merge_results_task(const depthn_result& left, const depthn_result& right, depthn_result& result, std::mutex& debug_cerr_mutex);
        void merge_results_task(const depthn_result& one, const depthn_result& two, const depthn_result& three, depthn_result& result, std::mutex& debug_cerr_mutex);
        //merge depth1 result with depthn result, special case (unique file in last batch)
        void merge_results_task(const depthn_result& left, const Depth1Result& right, depthn_result& result, std::mutex& debug_cerr_mutex);
        void merge_results_task(const depthn_result& left, const depthn_result& right, colors_to_minmer& final_result, std::mutex& debug_cerr_mutex);

        void read_file_task(const std::string& file, uint32_t doc_id, Depth1Result& result, const build::options_t& build_parameters, std::mutex& debug_cerr_mutex); 

        void init_batch(
            std::stack<Function>& task_stack,
            std::mutex& task_stack_mutex,
            std::deque<Depth1Result>& d1_results_storage,
            std::mutex& d1_results_storage_mutex,
            uint32_t start,
            uint32_t end,
            std::atomic<uint32_t>& running_task,
            std::condition_variable& cv,
            const build::options_t& build_parameters,
            std::mutex& debug_cerr_mutex);
        void run_batch(
            std::uint16_t batch_id,
            std::stack<Function>& task_stack,
            std::mutex& task_stack_mutex,
            std::deque<Depth1Result>& d1_results_storage,
            std::mutex& d1_results_storage_mutex,
            std::deque<depthn_result>& results_storage,
            std::atomic<uint32_t>& running_task,
            std::condition_variable& cv,
            const build::options_t& build_parameters,
            std::mutex& debug_cerr_mutex);
        void process_files(
            std::stack<Function>& task_stack,
            std::mutex& task_stack_mutex,
            std::deque<depthn_result>& results_storage,
            std::condition_variable& cv,
            std::atomic<uint32_t>& running_task,
            std::mutex& debug_cerr_mutex,
            const build::options_t& build_parameters);

        void run_batch_tree(
            uint16_t batch_id,
            uint32_t batch_size,
            std::stack<Function>& task_stack,
            std::mutex& task_stack_mutex,
            std::deque<depthn_result>& results_storage,
            std::mutex& results_storage_mutex,
            std::condition_variable& cv,
            std::atomic<uint32_t>& running_task,
            bool final_batch,
            colors_to_minmer& final_result,
            const build::options_t& build_parameters,
            std::mutex& debug_cerr_mutex); 
        void process_tree(
            std::stack<Function>& task_stack,
            std::mutex& task_stack_mutex,
            std::deque<depthn_result>& results_storage,
            std::mutex& results_storage_mutex,
            std::condition_variable& cv,
            std::atomic<uint32_t>& running_task,
            colors_to_minmer& final_result,
            std::mutex& debug_cerr_mutex,
            const build::options_t& build_parameters);

        //following methods are explicitly instantiated in src/psa/files
        //with colorsclasses being from hybrid.hpp and color mapper being pthash::compact_vector
        std::vector<color_t> full_dense_intersection(std::vector<typename ColorClasses::row_accessor>&& color_id_itrs) const noexcept;
        std::vector<color_t> full_mixed_intersection(std::vector<typename ColorClasses::row_accessor>&& color_id_itrs) const noexcept;

        std::vector<scored_id> ranking_dense_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold, uint64_t nb_kmers) const noexcept;
        std::vector<scored_id> ranking_mixed_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold) const noexcept;

        std::vector<color_t> union_dense_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold, uint64_t nb_kmers) const noexcept;
        std::vector<color_t> union_mixed_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold) const noexcept;
        
        std::vector<std::string> m_filenames;
        uint8_t k;
        uint8_t m;
        uint64_t seed;
        bool canonical;
        double pthash_constant;
        ColorClasses m_ccs; // colors
        ColorMapper m_map; // map between mphf values and color classes
        pthash_minimizers_mphf_t hf; // minimizer mphf
};

CLASS_HEADER
METHOD_HEADER::index()
    :
    k(0),
    m(0),
    seed(0),
    canonical(false),
    pthash_constant(0)
{}

CLASS_HEADER
METHOD_HEADER::index(const build::options_t& build_parameters)
    : 
    m_filenames(build_parameters.input_filenames),
    k(build_parameters.k),
    m(build_parameters.m),
    seed(build_parameters.seed),
    canonical(build_parameters.canonical),
    pthash_constant(build_parameters.pthash_constant)
{
    build2(build_parameters);
}



CLASS_HEADER
void 
METHOD_HEADER::memory_breakdown(std::ostream& out) const noexcept
{
    libra scale;
    scale.visit(m_filenames);
    out << "The list of input filenames weights: " << scale.get_byte_size() << " Bytes\n";
    scale.reset();
    out << "The MPHF of minimizers weights: " << hf.num_bits() / 8 << " Bytes\n";
    scale.visit(m_ccs);
    out << "Colors weight: " << scale.get_byte_size() << " Bytes\n";
    scale.reset();
    //TODO:fix this
    //scale.visit(m_map);
    out << "The mapping from minimizers to colors weights: " << m_map.bytes() << " Bytes\n";
    //scale.reset();
}

CLASS_HEADER
template <class Visitor>
void
METHOD_HEADER::visit(Visitor& visitor)
{
    visitor.visit(m_filenames);
    visitor.visit(k);
    visitor.visit(m);
    visitor.visit(seed);
    visitor.visit(canonical);
    visitor.visit(pthash_constant);
    visitor.visit(hf); // lphash mphf
    visitor.visit(m_ccs); // colors
    visitor.visit(m_map); // map between mphf values and color classes
}

CLASS_HEADER
typename METHOD_HEADER::pthash_opt_t
METHOD_HEADER::get_pthash_options(const build::options_t& build_parameters)
{
    pthash_opt_t opts;
    opts.seed = build_parameters.seed;
    opts.lambda = build_parameters.pthash_constant; // (too slow = try decreasing), higher lambda : more space efficient 
    opts.alpha = 0.97; //was 0.94
    opts.search = pthash::pthash_search_type::add_displacement;
    opts.avg_partition_size = 3000;
    opts.minimal_output = true;
    opts.verbose_output = build_parameters.verbose;
    
    opts.ram = build_parameters.max_ram * constants::GB;
    opts.num_threads = build_parameters.nthreads;
    opts.tmp_dir = build_parameters.tmp_dir;

    opts.dense_partitioning = true;

    return opts;
}



/*
Build 2 is a parallel version of build, using a thread pool to read files and merge results.
*/

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


// Worker thread function
CLASS_HEADER
void
METHOD_HEADER::worker_thread(
    std::stack<Function>& task_stack,
    std::mutex& task_stack_mutex,
    std::condition_variable& cv,
    std::atomic<uint32_t>& running_task,
    std::atomic<bool>& all_done,
    std::mutex& debug_cerr_mutex) 
{
    while (true) {
        Function task;
        {
            std::unique_lock<std::mutex> lock(task_stack_mutex);
            cv.wait(lock, [&]() {
                return !task_stack.empty() || all_done.load(); 
            });

            if (all_done.load() && task_stack.empty()) { //all done should be set to 1 only if stack empty but w/e
                break; // Exit if no tasks and manager is done
            }

            if (!task_stack.empty()) {
                running_task++; //already one in task but security (2 tasks running for each task)
                task = std::move(task_stack.top());
                task_stack.pop();
            }
        }
        if (task.f) {
            {
                std::lock_guard<std::mutex> lock(debug_cerr_mutex);
                std::cerr << "Running task: " << task.desc << std::endl;
            }
            task.f();
            running_task--;
            cv.notify_all();
        }
    }
}


CLASS_HEADER
void
METHOD_HEADER::init_batch(
    std::stack<Function>& task_stack,
    std::mutex& task_stack_mutex,
    std::deque<Depth1Result>& d1_results_storage,
    std::mutex& d1_results_storage_mutex,
    uint32_t start,
    uint32_t end,
    std::atomic<uint32_t>& running_task,
    std::condition_variable& cv,
    const build::options_t& build_parameters,
    std::mutex& debug_cerr_mutex)
{
    uint64_t allocated_ram = build_parameters.max_ram * constants::GB / 2 / (end - start);
    std::cerr << "Allocated RAM: " << allocated_ram << " Bytes\n";
    {
        std::lock_guard<std::mutex> lock(task_stack_mutex);
        for (uint32_t docid = start; docid < end; ++docid) {
            //create read_file_task for each file of batch
            task_stack.push( //Function(f, desc)
                {
                    [this, &start, &allocated_ram, docid, &d1_results_storage, &d1_results_storage_mutex, &running_task, &cv, &build_parameters, &debug_cerr_mutex]() {
                        ++running_task;
                        {
                            std::lock_guard<std::mutex> lock(debug_cerr_mutex);
                            std::cerr << "gonna create a d1 result with allocated ram: " << allocated_ram << " Bytes\n";
                        }
                        Depth1Result result(
                            {
                                emem_minmers(
                                    allocated_ram, 
                                    build_parameters.tmp_dir, 
                                    utils::get_tmp_filename("minmers", start, 0, std::this_thread::get_id())),
                                docid
                            }
                        );
                                               
                        read_file_task(m_filenames[docid], docid, result, build_parameters, debug_cerr_mutex);
                        {
                            std::lock_guard<std::mutex> lock(d1_results_storage_mutex);
                            d1_results_storage.push_back(std::move(result));
                        }

                        --running_task;
                        cv.notify_all();
                    },
                    "read_file_task" //function desc
                }
            );
        }
    } //end lock on stack
    cv.notify_all(); // Notify threads about new tasks
}


CLASS_HEADER
void
METHOD_HEADER::run_batch(
    uint16_t batch_id,
    std::stack<Function>& task_stack,
    std::mutex& task_stack_mutex,
    std::deque<Depth1Result>& d1_results_storage,
    std::mutex& d1_results_storage_mutex,
    std::deque<depthn_result>& results_storage,
    std::atomic<uint32_t>& running_task,
    std::condition_variable& cv,
    const build::options_t& build_parameters,
    std::mutex& debug_cerr_mutex)
{
    std::deque<depthn_result> tmp_results_storage;
    std::mutex tmp_results_storage_mutex;
    uint16_t depth = 1;
    uint64_t ram_for_1_depth = build_parameters.max_ram * constants::GB / 2;
    uint64_t allocated_ram;
    uint32_t nb_d1 = d1_results_storage.size();

    if (nb_d1 == 1 && m_filenames.size() == 1) {
        std::cerr << "TODO: handle case 1 file in index, likely going to crash\n";
    }
    //1ST STEP : PROCESS DEPTH 1 RESULTS
    allocated_ram = (nb_d1 == 1) ? ram_for_1_depth : ram_for_1_depth / (nb_d1/2); //divide by the nb of results in next depth

    if (nb_d1 % 2 == 1) {
        //odd number of files in batch, possibly last batch (if even number of threads)
        if (nb_d1 == 1 && m_filenames.size() > 1) {
            //case last batch contains only one file, merge it with last batch result
            depthn_result result(
                    ram_for_1_depth, 
                    build_parameters.tmp_dir, 
                    utils::get_tmp_filename("result", batch_id-1, -1, std::this_thread::get_id())
            );
            merge_results_task(results_storage.back(), d1_results_storage.front(), result, debug_cerr_mutex);
            result.minimize();
            results_storage.pop_back();
            results_storage.push_back(std::move(result));
            return;
        }

        else {
            //3 files in last batch or odd number (make it even by merging 3 firsts together)
            std::lock_guard<std::mutex> lock(task_stack_mutex);
            task_stack.push(
                {
                    [this, &batch_id, &depth, &allocated_ram, &d1_results_storage, &d1_results_storage_mutex, &tmp_results_storage, &tmp_results_storage_mutex, &build_parameters, &cv, &running_task, &debug_cerr_mutex]() {
                        //merge 3 results
                        ++running_task;
                        Depth1Result left, middle, right;
                        {
                            std::lock_guard<std::mutex> lock(d1_results_storage_mutex);
                            left = std::move(d1_results_storage.front());
                            d1_results_storage.pop_front();
                            middle = std::move(d1_results_storage.front());
                            d1_results_storage.pop_front();
                            right = std::move(d1_results_storage.front());
                            d1_results_storage.pop_front();
                        }
                    
                        depthn_result result(
                                allocated_ram, //divide by the nb of results in next depth
                                build_parameters.tmp_dir, 
                                utils::get_tmp_filename("result", batch_id, depth, std::this_thread::get_id())
                        );
                        merge_results_task(left, middle, right, result, debug_cerr_mutex);
                        {
                            std::lock_guard<std::mutex> lock(tmp_results_storage_mutex);
                            tmp_results_storage.push_back(std::move(result));
                        }
                        --running_task;
                        cv.notify_all();
                    },
                    "merge31_results_task"
                }
            );   
            nb_d1 = nb_d1 - 3; 
        }
    } //end case odd number of files in batch

    { //even number of files in batch, merge them 2 by 2
        std::lock_guard<std::mutex> lock(task_stack_mutex);
        while (nb_d1 != 0) {
            task_stack.push(
                {
                    [this, &batch_id, &depth, &allocated_ram, &d1_results_storage, &d1_results_storage_mutex, &tmp_results_storage, &tmp_results_storage_mutex, &build_parameters, &cv, &running_task, &debug_cerr_mutex]() {
                        ++running_task;
                        Depth1Result left, right;
                        {
                            std::lock_guard<std::mutex> lock(d1_results_storage_mutex);

                            left = std::move(d1_results_storage.front());
                            d1_results_storage.pop_front();

                            right = std::move(d1_results_storage.front());
                            d1_results_storage.pop_front();
                        }

                        // Create result object for merging
                        depthn_result result(
                            allocated_ram, 
                            build_parameters.tmp_dir, 
                            utils::get_tmp_filename("result", batch_id, depth, std::this_thread::get_id())
                        );

                        // Perform the merge
                        merge_results_task(left, right, result, debug_cerr_mutex);

                        {
                            std::lock_guard<std::mutex> lock(tmp_results_storage_mutex);
                            tmp_results_storage.push_back(std::move(result));
                        }
                        --running_task;
                        cv.notify_all();
                    },
                    "merge21_results_task"
                }
            );
            nb_d1 = nb_d1 - 2;
        }
    }

    //all merging tasks added to stack, notify threads
    cv.notify_all();

    //wait for all tasks to finish
    {
        std::unique_lock<std::mutex> lock(task_stack_mutex);
        cv.wait(lock, [&]() {
            return running_task.load() == 0 && task_stack.empty();
        });
    }
    
    //2ND STEP : MERGE DEPTH N RESULTS
    depth++;
    uint32_t nb_results = tmp_results_storage.size();
     //divide by the nb of results in next depth
    while (nb_results != 1){
        allocated_ram = (nb_results == 1) ? ram_for_1_depth : ram_for_1_depth / (nb_results/2);
        //iteratively go deeper into batch until everything merged 
        {
            std::lock_guard<std::mutex> lock(task_stack_mutex);
            if (nb_results % 2 == 1) {
                //odd number of results, 3 firsts are merged together
                task_stack.push(
                    {
                        [this, &batch_id, &depth, &allocated_ram, &tmp_results_storage, &tmp_results_storage_mutex, &build_parameters, &cv, &running_task, &debug_cerr_mutex]() {
                            ++running_task;
                            depthn_result left, middle, right;
                            {
                                std::lock_guard<std::mutex> lock(tmp_results_storage_mutex);
                                left = std::move(tmp_results_storage.front());
                                tmp_results_storage.pop_front();
                                middle = std::move(tmp_results_storage.front());
                                tmp_results_storage.pop_front();
                                right = std::move(tmp_results_storage.front());
                                tmp_results_storage.pop_front();
                            }

                            depthn_result result(
                                    allocated_ram, //divide by the nb of results in next depth
                                    build_parameters.tmp_dir, 
                                    utils::get_tmp_filename("result", batch_id, depth, std::this_thread::get_id())
                            );
                            merge_results_task(left, middle, right, result, debug_cerr_mutex);
                            {
                                std::lock_guard<std::mutex> lock(tmp_results_storage_mutex);
                                tmp_results_storage.push_back(std::move(result));
                            }
                            --running_task;
                            cv.notify_all();
                        },
                        "merge3n_results_task"
                    }
                );
                nb_results = nb_results - 3;
            }

            while (nb_results != 0) {
                //even number of results, merge them 2 by 2
                task_stack.push(
                    {
                        [this, &batch_id, &depth, &allocated_ram, &tmp_results_storage, &tmp_results_storage_mutex, &build_parameters, &cv, &running_task, &debug_cerr_mutex]() {
                            ++running_task;
                            depthn_result left, right;
                            {
                                std::lock_guard<std::mutex> lock(tmp_results_storage_mutex);
                                left = std::move(tmp_results_storage.front());
                                tmp_results_storage.pop_front();
                                right = std::move(tmp_results_storage.front());
                                tmp_results_storage.pop_front();
                            }

                            depthn_result result(
                                    allocated_ram,
                                    build_parameters.tmp_dir, 
                                    utils::get_tmp_filename("result", batch_id, depth, std::this_thread::get_id())
                            );
                            merge_results_task(left, right, result, debug_cerr_mutex);
                            {
                                std::lock_guard<std::mutex> lock(tmp_results_storage_mutex);
                                tmp_results_storage.push_back(std::move(result));
                            }
                            --running_task;
                            cv.notify_all();
                        },
                        "merge2n_results_task"
                    }
                );
                nb_results = nb_results - 2;
            }
        } //end lock on stack
        cv.notify_all(); // Notify threads about new tasks

        //wait for all tasks to finish
        {
            std::unique_lock<std::mutex> lock(task_stack_mutex);
            cv.wait(lock, [&]() {
                return running_task.load() == 0 && task_stack.empty();
            });
        }
            

        nb_results = tmp_results_storage.size();
        depth++;
    }

    //all results merged, transfer to results_storage
    tmp_results_storage.front().minimize();
    results_storage.push_back(std::move(tmp_results_storage.front()));
}


CLASS_HEADER
void
METHOD_HEADER::process_files(
    std::stack<Function>& task_stack,
    std::mutex& task_stack_mutex,
    std::deque<depthn_result>& results_storage,
    std::condition_variable& cv,
    std::atomic<uint32_t>& running_task,
    std::mutex& debug_cerr_mutex,
    const build::options_t& build_parameters)
{
    std::deque<Depth1Result> d1_results_storage;
    std::mutex d1_results_storage_mutex;
    uint16_t batch_id = 0;

    for (size_t i = 0; i < m_filenames.size(); i += build_parameters.nthreads) {
        //for each batch of files, do the parsing and merging
        uint32_t start = i;
        uint32_t end = (i + build_parameters.nthreads < m_filenames.size()) ? i + build_parameters.nthreads : m_filenames.size();

        std::cerr << "\n\n// BATCH " << batch_id << " //\n";

        init_batch(task_stack, task_stack_mutex, d1_results_storage, d1_results_storage_mutex, start, end, running_task, cv, build_parameters, debug_cerr_mutex);

        //wait for all tasks to finish
        {
            std::unique_lock<std::mutex> lock(task_stack_mutex);
            cv.wait(lock, [&]() {
                return running_task.load() == 0 && task_stack.empty();
            });
        }

        run_batch(batch_id, task_stack, task_stack_mutex, d1_results_storage, d1_results_storage_mutex, results_storage, running_task, cv, build_parameters, debug_cerr_mutex);

        batch_id++;
    }
}







CLASS_HEADER
void
METHOD_HEADER::run_batch_tree(
    uint16_t batch_id,
    uint32_t batch_size,
    std::stack<Function>& task_stack,
    std::mutex& task_stack_mutex,
    std::deque<depthn_result>& results_storage,
    std::mutex& results_storage_mutex,
    std::condition_variable& cv,
    std::atomic<uint32_t>& running_task,
    bool final_batch,
    colors_to_minmer& final_result,
    const build::options_t& build_parameters,
    std::mutex& debug_cerr_mutex)
{
    if (batch_size == 1){
        return;
    }

    std::deque<depthn_result> tmp_results_storage;
    std::mutex tmp_results_storage_mutex;
    uint16_t depth = build_parameters.nthreads;
    uint64_t ram_for_1_depth = build_parameters.max_ram * constants::GB / 2;
    uint64_t allocated_ram;

    uint32_t nb_results = batch_size;
    for (size_t i; i < batch_size; i++){
        tmp_results_storage.push_back(std::move(results_storage.front()));
        results_storage.pop_front();
    }

    if (final_batch){
        while (nb_results != 2){
            allocated_ram = (nb_results%2==0) ? (ram_for_1_depth / (nb_results/2)) : (ram_for_1_depth / ((nb_results+1)/2));
            {
                std::lock_guard<std::mutex> lock(task_stack_mutex);
                while (nb_results > 1) {
                    task_stack.push(
                        {
                            [this, &batch_id, &depth, &allocated_ram, &tmp_results_storage, &tmp_results_storage_mutex, &build_parameters, &cv, &running_task, &debug_cerr_mutex]() {
                                ++running_task;
                                depthn_result left, right;
                                {
                                    std::lock_guard<std::mutex> lock(tmp_results_storage_mutex);
                                    left = std::move(tmp_results_storage.front());
                                    tmp_results_storage.pop_front();
                                    right = std::move(tmp_results_storage.front());
                                    tmp_results_storage.pop_front();
                                }

                                depthn_result result(
                                        allocated_ram,
                                        build_parameters.tmp_dir, 
                                        utils::get_tmp_filename("tree_result", batch_id, depth, std::this_thread::get_id())
                                );
                                merge_results_task(left, right, result, debug_cerr_mutex);
                                {
                                    std::lock_guard<std::mutex> lock(tmp_results_storage_mutex);
                                    tmp_results_storage.push_back(std::move(result));
                                }
                                --running_task;
                                cv.notify_all();
                            },
                            "merge2n_final_task"
                        }
                    );
                    nb_results = nb_results - 2;
                }
            }
            cv.notify_all(); // Notify threads about new tasks

            //wait for all tasks to finish
            {
                std::unique_lock<std::mutex> lock(task_stack_mutex);
                cv.wait(lock, [&]() {
                    return running_task.load() == 0 && task_stack.empty();
                });
            }
            nb_results = tmp_results_storage.size();
            depth++;
        }
        // only 2 results left, merge them into a sorted emem, colors to minmers
        merge_results_task(tmp_results_storage.front(), tmp_results_storage.back(), final_result, debug_cerr_mutex);

        return;
    }

    //crash here if only one result in last batch, did special case for now
    while (nb_results != 1){
        allocated_ram = (nb_results == 1) ? ram_for_1_depth : ram_for_1_depth / (nb_results/2);
        //iteratively go deeper into batch until everything merged 
        {
            std::lock_guard<std::mutex> lock(task_stack_mutex);
            if (nb_results % 2 == 1) {
                //odd number of results, 3 firsts are merged together
                task_stack.push(
                    {
                        [this, &batch_id, &depth, &allocated_ram, &tmp_results_storage, &tmp_results_storage_mutex, &build_parameters, &cv, &running_task, &debug_cerr_mutex]() {
                            ++running_task;
                            depthn_result left, middle, right;
                            {
                                std::lock_guard<std::mutex> lock(tmp_results_storage_mutex);
                                left = std::move(tmp_results_storage.front());
                                tmp_results_storage.pop_front();
                                middle = std::move(tmp_results_storage.front());
                                tmp_results_storage.pop_front();
                                right = std::move(tmp_results_storage.front());
                                tmp_results_storage.pop_front();
                            }

                            depthn_result result(
                                    allocated_ram, //divide by the nb of results in next depth
                                    build_parameters.tmp_dir, 
                                    utils::get_tmp_filename("tree_result", batch_id, depth, std::this_thread::get_id())
                            );
                            merge_results_task(left, middle, right, result, debug_cerr_mutex);
                            {
                                std::lock_guard<std::mutex> lock(tmp_results_storage_mutex);
                                tmp_results_storage.push_back(std::move(result));
                            }
                            --running_task;
                            cv.notify_all();
                        },
                        "merge3n_results_task"
                    }
                );
                nb_results = nb_results - 3;
            }

            while (nb_results != 0) {
                //even number of results, merge them 2 by 2
                task_stack.push(
                    {
                        [this, &batch_id, &depth, &allocated_ram, &tmp_results_storage, &tmp_results_storage_mutex, &build_parameters, &cv, &running_task, &debug_cerr_mutex]() {
                            ++running_task;
                            depthn_result left, right;
                            {
                                std::lock_guard<std::mutex> lock(tmp_results_storage_mutex);
                                left = std::move(tmp_results_storage.front());
                                tmp_results_storage.pop_front();
                                right = std::move(tmp_results_storage.front());
                                tmp_results_storage.pop_front();
                            }

                            depthn_result result(
                                    allocated_ram,
                                    build_parameters.tmp_dir, 
                                    utils::get_tmp_filename("tree_result", batch_id, depth, std::this_thread::get_id())
                            );
                            merge_results_task(left, right, result, debug_cerr_mutex);
                            {
                                std::lock_guard<std::mutex> lock(tmp_results_storage_mutex);
                                tmp_results_storage.push_back(std::move(result));
                            }
                            --running_task;
                            cv.notify_all();
                        },
                        "merge2n_results_task"
                    }
                );
                nb_results = nb_results - 2;
            }
        } //end lock on stack
        cv.notify_all(); // Notify threads about new tasks

        //wait for all tasks to finish
        {
            std::unique_lock<std::mutex> lock(task_stack_mutex);
            cv.wait(lock, [&]() {
                return running_task.load() == 0 && task_stack.empty();
            });
        }
            

        nb_results = tmp_results_storage.size();
        depth++;
    }

    //all results merged, transfer to results_storage
    tmp_results_storage.front().minimize();
    results_storage.push_back(std::move(tmp_results_storage.front()));
}




CLASS_HEADER
void
METHOD_HEADER::process_tree(
    std::stack<Function>& task_stack,
    std::mutex& task_stack_mutex,
    std::deque<depthn_result>& results_storage,
    std::mutex& results_storage_mutex,
    std::condition_variable& cv,
    std::atomic<uint32_t>& running_task,
    colors_to_minmer& final_result,
    std::mutex& debug_cerr_mutex,
    const build::options_t& build_parameters)
{
    bool final_batch = false;
    uint16_t batch_id = 0;

    while (results_storage.size() != 1){
        uint32_t nb_results = results_storage.size();
        final_batch = (nb_results <= build_parameters.nthreads);

        for (size_t i = 0; i < nb_results; i += build_parameters.nthreads){
            std::cerr << "\n\n// BATCH " << batch_id << " //\n";
            //for each batch of results, do the merging
            uint32_t batch_size = (i + build_parameters.nthreads < nb_results) ? build_parameters.nthreads : nb_results-i;

            run_batch_tree(batch_id, batch_size, task_stack, task_stack_mutex,  results_storage, results_storage_mutex, cv, running_task, final_batch, final_result, build_parameters, debug_cerr_mutex);

            batch_id++;
            if (final_batch) return;
        }
    }
}




CLASS_HEADER
void
METHOD_HEADER::build2(const build::options_t& build_parameters)
{
    std::cerr << "DEBUG BUILD 2\n";
    auto start_time = std::chrono::high_resolution_clock::now();

    assert(m_filenames.size() > build_parameters.nthreads);
    //TODO, allow final_result in process_files(), for now only in process_tree()

    //synchro variables
    std::condition_variable cv;
    std::atomic<uint32_t> running_task(0);
    std::atomic<bool> all_done(false);

    std::stack<Function> task_stack;
    std::mutex task_stack_mutex;

    std::deque<depthn_result> results_storage;
    std::mutex results_storage_mutex;

    colors_to_minmer final_result(
        build_parameters.max_ram * constants::GB / 2, 
        build_parameters.tmp_dir, 
        utils::get_tmp_filename("final_result", 0, 0, std::this_thread::get_id())
    );

    std::mutex debug_cerr_mutex;

    // Create worker threads
    std::vector<std::thread> workers;
    for (uint32_t i = 0; i < build_parameters.nthreads; ++i) {
        workers.push_back(std::thread(
            &METHOD_HEADER::worker_thread, 
            this, std::ref(task_stack), std::ref(task_stack_mutex), std::ref(cv), std::ref(running_task), std::ref(all_done), std::ref(debug_cerr_mutex)
        ));
    }

    process_files(task_stack, task_stack_mutex, results_storage, cv, running_task, debug_cerr_mutex, build_parameters);

    {
        std::lock_guard<std::mutex> lock(debug_cerr_mutex);
        std::cerr << "All files processed\n";
        std::cerr << results_storage.size() << " results in storage\n";
        std::cerr << results_storage.front().size() << " results in first depthn\n";
    }

    std::cout << "DEBUG Time for process_files: " 
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - start_time
                    ).count() 
                << " milliseconds\n";

    start_time = std::chrono::high_resolution_clock::now();

    process_tree(task_stack, task_stack_mutex, results_storage, results_storage_mutex, cv, running_task, final_result, debug_cerr_mutex, build_parameters);

    // Mark all tasks as done
    all_done.store(true);
    cv.notify_all();

    {
        std::lock_guard<std::mutex> lock(debug_cerr_mutex);
        std::cerr << "All tree processed\n";
        std::cerr << results_storage.size() << " results in storage\n";
        std::cerr << final_result.size() << " pairs in final\n";
    }

    // Wait for worker threads to finish
    for (auto& worker : workers) {
        worker.join();
    }

    std::cout << "DEBUG Time for process_tree: " 
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - start_time
                    ).count() 
                << " milliseconds\n";


    ankerl::unordered_dense::set<uint64_t> unique_minmers;
    //TODO : can be done inside process_tree() to avoid running through again
    auto itr = final_result.cbegin();
    while(itr != final_result.cend()) {
        unique_minmers.insert((*itr).second);
        ++itr;
    }

    {
    //STEP 3 : BUILDING MPHF ===================================================
        if (build_parameters.verbose > 0) std::cerr << "Step 3: building the MPHF for " << unique_minmers.size() << " minimizers\n";
        start_time = std::chrono::high_resolution_clock::now();
        //auto pt_itr = pthash_input_iterator<decltype(unique_minimizers)::const_iterator>(unique_minimizers.cbegin());
        int backup, redirect;
        fflush(stdout);
        backup = dup(1);
        redirect = open("/dev/null", O_WRONLY);
        dup2(redirect, 1);
        close(redirect);

        //MPHF via PTHash, n minmers, hf : minmer(uint64) -> [0, n-1]
        hf.build_in_internal_memory(unique_minmers.begin(), unique_minmers.size(), get_pthash_options(build_parameters));

        fflush(stdout);
        dup2(backup, 1);
        close(backup);
        assert(hf.num_keys() == unique_minmers.size());
        unique_minmers.clear();

        std::cout << "DEBUG Time for MPHF build: " 
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - start_time
                    ).count() 
                << " milliseconds\n";

        std::cerr << "DEBUG End building MPHF\n";
        std::cerr << "size of MPHF : " << hf.num_bits()/8 << "\n";
        utils::printRAMInfo();
    }

    std::cerr << "DEBUG in-between step 3 and 4\n";
    utils::printRAMInfo();

    //STEP 4 : LIST DEDUPLICATION + MAPPING ====================================
    {
        if (build_parameters.verbose > 0) std::cerr << "Step 4: list deduplication and mapping\n";
        start_time = std::chrono::high_resolution_clock::now();
        
        typename ColorClasses::builder cbuild(m_filenames.size(), build_parameters.verbose);
        pthash::compact_vector::builder m_map_builder(hf.num_keys(), ceil(log2(hf.num_keys())));

        if (build_parameters.check) {
            if (build_parameters.verbose > 0) std::cerr << "map/MPHF of size: " << hf.num_keys() << "\n";
            for (std::size_t i = 0; i < hf.num_keys(); ++i) 
                m_map_builder.set(i, (1 << m_map_builder.width()) - 1); //max value with ceil(log2(hf.num_keys())) bits
        }

        uint32_t cid = 0;
        auto itr = final_result.cbegin();
        while(itr != final_result.cend()) {
            auto current_color = (*itr).first;
            cbuild.add_color_set(current_color.data(), current_color.size()); // only one copy of the color in storage
            while(itr != final_result.cend() and (*itr).first == current_color) {
                auto minimizer = (*itr).second;
                auto mp_idx = hf(minimizer);
                // std::cerr << minimizer << " -> " << mp_idx << "\n";
                m_map_builder.set(mp_idx, cid);
                ++itr;
            }
            ++cid;
        }
        cbuild.build(m_ccs);
        m_map_builder.build(m_map);
        //compact_vector m_map where m_map[ hf(minmer) ] = color_id 
    }

    std::cout << "DEBUG Time for dedup + mapping: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::high_resolution_clock::now() - start_time
                 ).count() 
              << " milliseconds\n";

    std::cerr << "DEBUG end mapping\n";
    std::cerr << "check mem breakdown for thr rest \n";
    utils::printRAMInfo();

    if (build_parameters.verbose > 0) {
        std::cerr << "Number of ids: " << m_ccs.num_docs() << "\n";
        std::cerr << "Number of colors (lists of ids):" << m_ccs.num_color_classes() << "\n";
    } 

}


#undef CLASS_HEADER
#undef METHOD_HEADER

} // namespace minimizer
} // namespace kaminari

#endif // KAMINARI_INDEX_HPP



/*Number of ids: 50
Number of colors (lists of ids):188898
The list of input filenames weights: 4134 Bytes
The MPHF of minimizers weights: 1059855 Bytes
Colors weight: 1708436 Bytes
The mapping from minimizers to colors weights: 8883056 Bytes

Written 11655515 Bytes
*/