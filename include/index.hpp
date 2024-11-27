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
typedef ankerl::unordered_dense::map<uint64_t, std::vector<uint32_t>> result_map;
typedef std::pair<ankerl::unordered_dense::set<uint64_t>, uint32_t> minmer_set_docid;

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
        void manager_thread(
            std::stack<Function>& task_stack,
            std::mutex& task_stack_mutex,
            std::deque<minmer_set_docid>& results_depth1,
            std::mutex& results_depth1_mutex,
            std::deque<result_map>& results,
            std::mutex& results_mutex,
            std::condition_variable& cv,
            std::atomic<uint32_t>& running_task,
            std::atomic<bool>& all_done,
            std::mutex& debug_cerr_mutex);
        result_map merge_results(const result_map& left, const result_map& right, std::mutex& debug_cerr_mutex);
        result_map merge_results(const minmer_set_docid& left, const minmer_set_docid& right, std::mutex& debug_cerr_mutex);
        // Function to read a file
        void read_file_task(const std::string& file, uint32_t doc_id, minmer_set_docid& result, const build::options_t& build_parameters, std::mutex& debug_cerr_mutex); 

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


CLASS_HEADER
void
METHOD_HEADER::build(const build::options_t& build_parameters)
{
    std::size_t run_id = std::chrono::high_resolution_clock::now().time_since_epoch().count();

    typedef std::pair<minimizer_t, color_t> mpc_t;
    emem::external_memory_vector<mpc_t> tmp_sorted_storage(
        build_parameters.max_ram * constants::GB, 
        build_parameters.tmp_dir, 
        utils::get_tmp_filename("", "minimizer_unitig_id", run_id)
    );

    std::unordered_map<uint64_t, std::vector<color_t>> kmer_color_check;


    auto start_time = std::chrono::high_resolution_clock::now();

    //STEP 1 : READING FILES ===================================================
    {
        if (build_parameters.verbose > 0) std::cerr << "Step 1: reading files\n";

        std::cerr << "DEBUG Start reading files\n";
        utils::printRAMInfo();

        float fraction = 0.1;
        color_t id = 0;
        std::size_t total_kmers = 0;
        std::size_t total_mmers = 0;
        std::size_t total_minimizers = 0;
        gzFile fp = nullptr;
        kseq_t* seq = nullptr;
        std::vector<::minimizer::record_t> mms_buffer;
        for (auto filename : m_filenames) {
            if ((fp = gzopen(filename.c_str(), "r")) == NULL)
                throw std::runtime_error("Unable to open input file " + filename);
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
                    std::cerr << "\tRead contig:\n" 
                              << "\t\t" << contig_kmer_count << " k-mers\n"
                              << "\t\t" << contig_mmer_count << " m-mers\n"
                              << "\t\t" << mms_buffer.size() << " minimizers\n";
                }
                // remove duplicates and output to tmp file on disk
                std::vector<minimizer_t> minimizers;
                for (auto r : mms_buffer) {
                    minimizers.push_back(r.itself);
                }
                std::sort(minimizers.begin(), minimizers.end());
                auto last = std::unique(minimizers.begin(), minimizers.end());
                for (auto itr = minimizers.begin(); itr != last; ++itr) {
                    tmp_sorted_storage.push_back(std::make_pair(*itr, id));
                }
                mms_buffer.clear();

                if (build_parameters.check) {
                    std::vector<::minimizer::record_t> kmer_buffer;
                    std::size_t dummy;
                    [[maybe_unused]] auto kmc = ::minimizer::from_string<hash64>(
                        seq->seq.s,
                        seq->seq.l,
                        build_parameters.k,
                        build_parameters.k,
                        static_cast<uint64_t>(0),
                        false, // minimizers are canonical or not, kmers are kept as they are for checking
                        dummy,
                        kmer_buffer
                    );
                    // auto last = std::unique(kmer_buffer.begin(), kmer_buffer.end());
                    for (auto itr = kmer_buffer.begin(); itr != kmer_buffer.end(); ++itr) {
                        if (kmer_color_check.find(itr->itself) == kmer_color_check.end()) {
                            kmer_color_check.emplace(itr->itself, std::vector<color_t> {id});
                        } else if (kmer_color_check[itr->itself].back() != id) { // this takes care of duplicates
                            kmer_color_check[itr->itself].push_back(id);
                        }
                    }
                }
            }
            if (seq) kseq_destroy(seq);
            gzclose(fp);
            fp = nullptr;
            seq = nullptr;
            ++id;
            /**/
            if (build_parameters.verbose > 1) {
                if (id >= m_filenames.size() * fraction) {
                    std::cerr << "\tProcessed " << fraction * 100 << "\% of input files\n";
                    std::cerr << "\t\t(minimizer, color) pairs: " << tmp_sorted_storage.size() << "\n";
                    fraction += 0.1;
                }
            }
        }

        std::cout << "DEBUG Time for reading files: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::high_resolution_clock::now() - start_time
                 ).count() 
              << " milliseconds\n";

        std::cerr << "DEBUG End reading files\n";
        std::cerr << "size of emem1 : " << tmp_sorted_storage.size()*sizeof(mpc_t) << "\n";
        utils::printRAMInfo();

        if (build_parameters.verbose > 0) {
            std::cerr << "\ttotal k-mers: " << total_kmers << "\n";
            std::cerr << "\ttotal m-mers: " << total_mmers << "\n";
            std::cerr << "\ttotal minimizers (with repetitions): " << total_minimizers << "\n";
        }
    }

    typedef std::pair<std::vector<color_t>, minimizer_t> cc_mm_t;
    emem::external_memory_vector<cc_mm_t> sorted_color_lists(
        build_parameters.max_ram * constants::GB, 
        build_parameters.tmp_dir, 
        utils::get_tmp_filename("", "color_minimizer_list", run_id)
    );

    std::cerr << "DEBUG in-between step and 2, defined(and reserved buffer) emem2\n";
    utils::printRAMInfo();


    //STEP 2 : AGGREGATING COLORS ==============================================
    { 
        // aggregate colors into lists for each minimizer (and build the MPHF while doing so)
        if (build_parameters.verbose > 0) std::cerr << "Step 2: aggregating colors\n";
        start_time = std::chrono::high_resolution_clock::now();
        /*emem::external_memory_vector<minimizer_t, false> unique_minimizers(
            build_parameters.max_ram * constants::GB, 
            build_parameters.tmp_dir, 
            utils::get_tmp_filename("", "unique_minimizers", run_id)
        );*/
        std::vector<minimizer_t> unique_minimizers; //TODO pthash phobic dense_partitonned needs random access to the keys for parallelism. emem random access ? for now vector but might overflow RAM (347M minimisers for 60human docs -> 1.39GB RAM)
        
        auto itr = tmp_sorted_storage.cbegin();
        while (itr != tmp_sorted_storage.cend()) { // build list of colors for each minimizer
            auto current = (*itr).first;
            std::vector<color_t> ids;
            while (itr != tmp_sorted_storage.cend() and (*itr).first == current) {
                ids.push_back((*itr).second);
                ++itr;
            }
            auto last = std::unique(ids.begin(), ids.end()); //do we expect any duplicates here? -> no because we already deduplicated the minimizers for each file, TODO: choice between deduplicating minimisers in one file here or before
            ids.erase(last, ids.end());
            sorted_color_lists.push_back(std::make_pair(std::move(ids), current)); // save ([ids], minimizer) to disk
            unique_minimizers.push_back(current);
        }

        std::cout << "DEBUG Time for aggregating colors: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::high_resolution_clock::now() - start_time
                 ).count() 
              << " milliseconds\n";

        std::cerr << "DEBUG End aggregating colors\n";
        std::cerr << "size of emem1 (should be deleted now) : " << tmp_sorted_storage.size()*sizeof(mpc_t) << "\n";
        std::cerr << "size of emem2 : " << sorted_color_lists.size()*sizeof(cc_mm_t) << "\n";
        std::cerr << "size of unique_minmers : " << unique_minimizers.size()*sizeof(minimizer_t) << "\n";
        utils::printRAMInfo();


    //STEP 3 : BUILDING MPHF ===================================================
        if (build_parameters.verbose > 0) std::cerr << "Step 3: building the MPHF for " << unique_minimizers.size() << " minimizers\n";
        start_time = std::chrono::high_resolution_clock::now();
        //auto pt_itr = pthash_input_iterator<decltype(unique_minimizers)::const_iterator>(unique_minimizers.cbegin());
        int backup, redirect;
        fflush(stdout);
        backup = dup(1);
        redirect = open("/dev/null", O_WRONLY);
        dup2(redirect, 1);
        close(redirect);

        //MPHF via PTHash, n minmers, hf : minmer(uint64) -> [0, n-1]
        hf.build_in_internal_memory(unique_minimizers.begin(), unique_minimizers.size(), get_pthash_options(build_parameters));

        fflush(stdout);
        dup2(backup, 1);
        close(backup);
        assert(hf.num_keys() == unique_minimizers.size());

        std::cout << "DEBUG Time for MPHF build: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::high_resolution_clock::now() - start_time
                 ).count() 
              << " milliseconds\n";

        std::cerr << "DEBUG End building MPHF\n";
        std::cerr << "size of emem1 (should be deleted already) : " << tmp_sorted_storage.size()*sizeof(mpc_t) << "\n";
        std::cerr << "size of emem2 : " << sorted_color_lists.size()*sizeof(cc_mm_t) << "\n";
        std::cerr << "size of unique_minmers : " << unique_minimizers.size()*sizeof(minimizer_t) << "\n";
        std::cerr << "size of MPHF : " << hf.num_bits()/8 << "\n";
        utils::printRAMInfo();
    }

    std::cerr << "DEBUG in-between step 3 and 4\n";
    std::cerr << "size of emem1 (should be deleted already) : " << tmp_sorted_storage.size()*sizeof(mpc_t) << "\n";
    std::cerr << "size of emem2 : " << sorted_color_lists.size()*sizeof(cc_mm_t) << "\n";
    std::cerr << "size of MPHF : " << hf.num_bits()/8 << "\n";
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
        auto itr = sorted_color_lists.cbegin();
        while(itr != sorted_color_lists.cend()) {
            auto current_color = (*itr).first;
            cbuild.add_color_set(current_color.data(), current_color.size()); // only one copy of the color in storage
            while(itr != sorted_color_lists.cend() and (*itr).first == current_color) {
                auto minimizer = (*itr).second;
                auto mp_idx = hf(minimizer);
                // std::cerr << minimizer << " -> " << mp_idx << "\n";
                if (build_parameters.check) { //super super slow but how to access in builder ?
                    pthash::compact_vector check_m_map;
                    m_map_builder.build(check_m_map);
                    if (check_m_map[mp_idx] != static_cast<uint64_t>((1 << m_map_builder.width()) - 1)) throw std::runtime_error("[check fail] reassigning id of unique minimizer (the minimizer is not unique)");
                }
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
    std::cerr << "size of emem1 (should be deleted already) : " << tmp_sorted_storage.size()*sizeof(mpc_t) << "\n";
    std::cerr << "size of emem2 : " << sorted_color_lists.size()*sizeof(cc_mm_t) << "\n";
    std::cerr << "size of MPHF : " << hf.num_bits()/8 << "\n";
    std::cerr << "check mem breakdown for thr rest \n";
    utils::printRAMInfo();

    if (build_parameters.verbose > 0) {
        std::cerr << "Number of ids: " << m_ccs.num_docs() << "\n";
        std::cerr << "Number of colors (lists of ids):" << m_ccs.num_color_classes() << "\n";
    }
}

/*
Build 2 is a parallel version of build, using a thread pool to read files and merge results.
*/

// Function to read a file
CLASS_HEADER
void 
METHOD_HEADER::read_file_task(const std::string& file, uint32_t doc_id, minmer_set_docid& result, const build::options_t& build_parameters, std::mutex& debug_cerr_mutex) {

    {
        std::lock_guard<std::mutex> lock(debug_cerr_mutex);
        std::cerr << "DEBUG Reading file " << file << "\n";
    }

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
            std::cerr << "\tRead contig:\n" 
                        << "\t\t" << contig_kmer_count << " k-mers\n"
                        << "\t\t" << contig_mmer_count << " m-mers\n"
                        << "\t\t" << mms_buffer.size() << " minimizers\n";
        }
        // remove duplicates and output to tmp file on disk
        std::vector<minimizer_t> minimizers;
        for (auto r : mms_buffer) {
            result.first.insert(r.itself);
        }
        mms_buffer.clear();
    }
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    fp = nullptr;
    seq = nullptr;
}

// Function to merge two depth_1 results
CLASS_HEADER
result_map
METHOD_HEADER::merge_results(const minmer_set_docid& left, const minmer_set_docid& right, std::mutex& debug_cerr_mutex) {
    result_map merged;
    // Determine insertion order based on docid
    if (left.second < right.second) {
        for (const auto& key : left.first) merged[key].push_back(left.second);
        for (const auto& key : right.first) merged[key].push_back(right.second);
    } else {
        for (const auto& key : right.first) merged[key].push_back(right.second);
        for (const auto& key : left.first) merged[key].push_back(left.second);
    }
    return merged;
}
// Function to merge two results
CLASS_HEADER
result_map
METHOD_HEADER::merge_results(const result_map& left, const result_map& right, std::mutex& debug_cerr_mutex) {
    result_map merged;
    for (const auto& [key, value] : left) {
        auto& vec = merged[key];
        vec.insert(vec.end(), value.begin(), value.end());
    }
    for (const auto& [key, value] : right) {
        auto& vec = merged[key];
        vec.insert(vec.end(), value.begin(), value.end());
    }
    for (auto& [key, vec] : merged) {
        std::sort(vec.begin(), vec.end()); //could be replaced by a set ? and make it a vector sorted at the end
        //vec.erase(std::unique(vec.begin(), vec.end()), vec.end()); -> useless for now because 1file : 1result, might change with future optis
    }
    return merged;
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
                task = std::move(task_stack.top());
                task_stack.pop();
            }
        }
        if (task.f) {
            task.f();
        }
    }
}


// Manager thread function
CLASS_HEADER
void
METHOD_HEADER::manager_thread(
    std::stack<Function>& task_stack,
    std::mutex& task_stack_mutex,
    std::deque<minmer_set_docid>& results_depth1,
    std::mutex& results_depth1_mutex,
    std::deque<result_map>& results,
    std::mutex& results_mutex,
    std::condition_variable& cv,
    std::atomic<uint32_t>& running_task,
    std::atomic<bool>& all_done,
    std::mutex& debug_cerr_mutex)
{
    while (true) {
        {
            std::lock_guard<std::mutex> resd1_lock(results_depth1_mutex);
            std::lock_guard<std::mutex> res_lock(results_mutex);
            std::lock_guard<std::mutex> stack_lock(task_stack_mutex);
            if (task_stack.empty() && running_task.load() == 0 && results_depth1.size() == 0 && results.size() == 1) { 
                //need no more tasks to do, no running tasks, 1 result left to get out
                all_done.store(true);
                cv.notify_all();

                // Double-check the state after setting all_done
                std::this_thread::sleep_for(std::chrono::milliseconds(250));
                if (!task_stack.empty() || running_task.load() > 0) {
                    all_done.store(false); // Revert if inconsistent
                } else {
                    break; // Consistent, safe to exit
                }
            }
        } //end of locks
    
        //main : look for merge to do
        {
            std::lock_guard<std::mutex> lock(results_mutex);
            // Find two results with the same depth to merge
            while (results.size() > 1) {
                {
                    std::lock_guard<std::mutex> lock(debug_cerr_mutex);
                    std::cerr << "DEBUG Merging results\n";
                }
                result_map left = std::move(results.front());
                results.pop_front();

                result_map right = std::move(results.front());
                results.pop_front();

                // Add merge task to the stack
                {
                    std::lock_guard<std::mutex> lock(task_stack_mutex);
                    
                    task_stack.push( //Function(f, desc)
                        {
                            [this, &running_task, left = std::move(left), right = std::move(right), &results, &results_mutex, &cv, &debug_cerr_mutex]() mutable {
                                ++running_task;
                                result_map merged = this->merge_results(left, right, debug_cerr_mutex);
                                {
                                    std::lock_guard<std::mutex> lock(results_mutex);
                                    results.push_back(std::move(merged));
                                }
                                --running_task;
                                cv.notify_all(); // Notify manager about new result
                            },
                            "merge_results_task" //function desc
                        }
                    );
                }
                cv.notify_all(); // Notify threads about new task
                //return;
            }
            //do the merging of depth_1 results after because want them on top of stack
            {
                std::lock_guard<std::mutex> lock(results_depth1_mutex);
                while (results_depth1.size() > 1) {
                    {
                        std::lock_guard<std::mutex> lock(debug_cerr_mutex);
                        std::cerr << "DEBUG Merging d1 results queue size : " << results_depth1.size() << "\n";
                    }
                    minmer_set_docid left = std::move(results_depth1.front());
                    results_depth1.pop_front();

                    minmer_set_docid right = std::move(results_depth1.front());
                    results_depth1.pop_front();

                    // Add merge task to the stack
                    {
                        std::lock_guard<std::mutex> lock(task_stack_mutex);
                        
                        task_stack.push( //Function(f, desc)
                            {
                                [this, &running_task, left = std::move(left), right = std::move(right), &results, &results_mutex, &cv, &debug_cerr_mutex]() mutable {
                                    ++running_task;
                                    result_map merged = this->merge_results(left, right, debug_cerr_mutex);
                                    {
                                        std::lock_guard<std::mutex> lock(results_mutex);
                                        results.push_back(std::move(merged));
                                    }
                                    --running_task;
                                    cv.notify_all(); // Notify manager about new result
                                },
                                "merge_results_d1_task" //function desc
                            }
                        );
                    }
                    cv.notify_all(); // Notify threads about new task
                } //end while
            } //end lock results_depth1_mutex
        } //end lock results_mutex 

        // TODO still miss odd number of files case : one d1 result left, need to merge it with a map


        //end : sleep until new results

        //TODO : dunno what to wait for, just wait some ms for now
        /* {
            std::unique_lock<std::mutex> lock(results_mutex);
            cv.wait(lock, [&]() { 
                return results.size() > 0;
            });
        } */
        //wait a bit to do not spent time looping and keep locking access to stack and results
        std::this_thread::sleep_for(std::chrono::milliseconds(5)); 
    }
}

CLASS_HEADER
void
METHOD_HEADER::build2(const build::options_t& build_parameters)
{
    std::cerr << "DEBUG BUILD 2\n";

    auto start_time = std::chrono::high_resolution_clock::now();

    std::mutex debug_cerr_mutex;

    // Shared state
    std::stack<Function> task_stack;
    std::mutex task_stack_mutex;

    std::deque<minmer_set_docid> results_depth1;
    std::mutex results_depth1_mutex;
    std::deque<result_map> results;
    std::mutex results_mutex;

    std::condition_variable cv;
    std::atomic<uint32_t> running_task(0);
    std::atomic<bool> all_done(false);

    // Initialize task stack with file-reading tasks
    int doc_id = 0;
    for (const auto& file : m_filenames) {
        task_stack.push( //Function(f, desc)
            {
                [this, &running_task, file, doc_id, &results_depth1, &results_depth1_mutex, &cv, &build_parameters, &debug_cerr_mutex]() mutable {
                    ++running_task;
                    minmer_set_docid result = {ankerl::unordered_dense::set<minimizer_t>(), doc_id};
                    this->read_file_task(file, doc_id, result, build_parameters, debug_cerr_mutex);
                    {
                        std::lock_guard<std::mutex> lock(results_depth1_mutex);
                        results_depth1.push_back(std::move(result));
                    }
                    --running_task;
                    cv.notify_all(); // Notify manager that a new result is available
                }, 
                "read_file_task" //function desc 
            }
        );
        doc_id++;
    }

    std::cerr << "DEBUG in-between step 1 and 2\n";

    // Start manager thread
    std::thread manager([this, &task_stack, &task_stack_mutex, &results_depth1, &results_depth1_mutex, &results, &results_mutex, &cv, &running_task, &all_done, &debug_cerr_mutex]() {
        this->manager_thread(task_stack, task_stack_mutex, results_depth1, results_depth1_mutex, results, results_mutex,  cv, running_task, all_done, debug_cerr_mutex);
    });

    std::cerr << "DEBUG gonna start my bois\n";

    // Start worker threads
    std::vector<std::thread> workers;
    for (int i = 0; i < build_parameters.nthreads - 1; ++i) {
        workers.emplace_back([this, &task_stack, &task_stack_mutex, &cv, &running_task, &all_done, &debug_cerr_mutex]() {
            this->worker_thread(task_stack, task_stack_mutex, cv, running_task, all_done,  debug_cerr_mutex);
        });
    }

    
    // Wait for workers to finish
    for (auto& worker : workers) {
        worker.join();
    }
    // Wait for manager to finish
    manager.join();

    assert(results.size() == 1);
    result_map& minmer_to_colors = results.front();
    ankerl::unordered_dense::set<uint64_t> unique_minmers;
    for (const auto& [minmer, colors] : minmer_to_colors) {
        unique_minmers.insert(minmer);
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

        std::unordered_map<std::vector<uint32_t>, uint32_t, VectorHash> color_to_cid;        
        uint32_t cid = 0;

        for (const auto& pair : minmer_to_colors) {
            const auto& minmer = pair.first;
            const auto& colors = pair.second;

            // Check if this vector already has a cid
            if (color_to_cid.find(colors) == color_to_cid.end()) {
                // Assign a new cid
                color_to_cid[colors] = cid++;
                cbuild.add_color_set(colors.data(), colors.size());
            }

            
            // Map minmer to the cid
            m_map_builder.set(hf(minmer), cid);
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
