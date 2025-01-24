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
#include "../bundled/biolib/include/iterator/standalone_iterator.hpp"
#include "../bundled/biolib/include/iterator/sorted_merge_iterator.hpp"
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

typedef std::pair<uint64_t, uint32_t> value_t;
typedef emem::external_memory_vector<value_t> emem_t;
typedef emem::external_memory_vector<std::pair<std::vector<uint32_t>, uint64_t>> colors_to_minmer;

struct Function {
    std::function<void()> f; // Actual function
    std::string desc;
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

        void worker_thread(
            std::stack<Function>& task_stack,
            std::mutex& task_stack_mutex,
            std::condition_variable& cv,
            std::atomic<uint32_t>& running_task,
            std::atomic<bool>& all_done);
        void read_file_task(const std::string& file, uint32_t doc_id, emem_t& result, const build::options_t& build_parameters); 
        void build(const build::options_t& build_parameters);


        //following methods are explicitly instantiated in src/psa/files
        //with colorsclasses being from hybrid.hpp and color mapper being pthash::compact_vector
        std::vector<scored_id> ranking_dense_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold) const noexcept;
        std::vector<scored_id> ranking_mixed_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold) const noexcept;

        std::vector<color_t> union_dense_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold) const noexcept;
        std::vector<color_t> union_mixed_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold) const noexcept;
        
        std::vector<std::string> m_filenames;
        uint8_t k;
        uint8_t m;
        uint8_t b;
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
    b(0),
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
    b(build_parameters.b),
    seed(build_parameters.seed),
    canonical(build_parameters.canonical),
    pthash_constant(build_parameters.pthash_constant)
{
    build(build_parameters);
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


// Worker thread function
CLASS_HEADER
void
METHOD_HEADER::worker_thread(
    std::stack<Function>& task_stack,
    std::mutex& task_stack_mutex,
    std::condition_variable& cv,
    std::atomic<uint32_t>& running_task,
    std::atomic<bool>& all_done) 
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
            task.f();
            running_task--;
            cv.notify_all();
        }
    }
}

// Function to read a file
CLASS_HEADER
void 
METHOD_HEADER::read_file_task(const std::string& file, uint32_t doc_id, emem_t& result, const build::options_t& build_parameters) {
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

        // duplicates removed during merge
        std::vector<minimizer_t> minimizers;
        for (auto r : mms_buffer) {
            result.push_back(
                std::make_pair(r.itself, doc_id));
        }
        mms_buffer.clear();
    }
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    fp = nullptr;
    seq = nullptr;
}

CLASS_HEADER
void
METHOD_HEADER::build(const build::options_t& build_parameters)
{
    //STEP 1 : PARSE FILES =====================================================
    if (build_parameters.verbose > 0) std::cerr << "Step 1: Parsing " << m_filenames.size() << " files\n";
    auto start_time = std::chrono::high_resolution_clock::now();
    
    //synchro variables
    std::condition_variable cv;
    std::atomic<uint32_t> running_task(0);
    std::atomic<bool> all_done(false);

    std::stack<Function> task_stack;
    std::mutex task_stack_mutex;

    std::vector<emem_t> results_storage;
    std::mutex results_mutex;

    
    //adding pase_file task to a stack of tasks to be picked by threads
    uint32_t doc_id = 0;
    for (const auto& file : m_filenames) {
        task_stack.push( //Function(f, desc)
            {
                [this, &running_task, file, doc_id, &results_storage, &results_mutex, &cv, &build_parameters]() mutable {
                    ++running_task;
                    emem_t result(
                        build_parameters.max_ram * constants::GB / build_parameters.nthreads,
                        build_parameters.tmp_dir, 
                        utils::get_tmp_filename("", "parse_result_doc", doc_id));

                    this->read_file_task(file, doc_id, result, build_parameters);

                    result.minimize();
                    {
                        std::lock_guard<std::mutex> lock(results_mutex);
                        results_storage.push_back(std::move(result));
                    }
                    --running_task;
                    cv.notify_all(); // Notify manager that a new result is available
                }, 
                "read_file_task" //function desc 
            }
        );
        doc_id++;
    }

    
    // Create worker threads who will pick tasks if there are any
    std::vector<std::thread> workers;
    for (uint32_t i = 0; i < build_parameters.nthreads; ++i) {
        workers.push_back(std::thread(
            &METHOD_HEADER::worker_thread, 
            this, std::ref(task_stack), std::ref(task_stack_mutex), std::ref(cv), std::ref(running_task), std::ref(all_done))
        );
    }

    //wait for all tasks to finish
    {
        std::unique_lock<std::mutex> lock(task_stack_mutex);
        cv.wait(lock, [&]() {
            return running_task.load() == 0 && task_stack.empty();
        });
    }

    // Mark all tasks as done
    all_done.store(true);
    cv.notify_all();

    for (auto& worker : workers) {
        worker.join();
    }

    std::cerr << "DEBUG Time for parsing: " 
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - start_time
                    ).count() 
                << " milliseconds\n";
    // all emem created, filled, sorted, dumped to disk

    //STEP 2 : MERGE RESULTS ===================================================
    if (build_parameters.verbose > 0) std::cerr << "Step 2: Merging results of parsing\n";
    start_time = std::chrono::high_resolution_clock::now();
    ankerl::unordered_dense::set<uint64_t> unique_minmers;

    colors_to_minmer final_result(
        build_parameters.max_ram * constants::GB * 0.75, //leave space for unique_minmers
        build_parameters.tmp_dir, 
        utils::get_tmp_filename("", "final_merge", 0)
    );


    std::vector<iterators::standalone_iterator<emem_t::const_iterator>> itr_vec;
    for (auto& emem : results_storage) {
        auto sa_itr = iterators::standalone::const_from(emem);
        itr_vec.push_back(sa_itr);
    }


    auto itr = iterators::sorted_merge_iterator(std::move(itr_vec));
    iterators::sorted_merge_iterator<emem_t::const_iterator> end;

    uint64_t smallest;
    while (itr != end){
        smallest = (*itr).first;
        std::pair<std::vector<uint32_t>, uint64_t> color_minmer = {{(*itr).second}, smallest};
        ++itr;
        while (itr != end && (*itr).first == smallest) {
            if ((*itr).second != color_minmer.first.back()){
                //remove duplicates inside 1 file
                color_minmer.first.push_back((*itr).second);
            }
            ++itr;
        }

        final_result.push_back(color_minmer);
        unique_minmers.insert(smallest);
    }

    final_result.minimize();

    std::cerr << "DEBUG Time for merging: " 
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - start_time
                    ).count() 
                << " milliseconds\n";

    //STEP 3 : BUILDING MPHF ===================================================
    {
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

        std::cerr << "DEBUG Time for MPHF build: " 
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - start_time
                    ).count() 
                << " milliseconds\n";

        std::cerr << "size of MPHF : " << hf.num_bits()/8 << "\n";
    }

    std::cerr << "DEBUG in-between step 3 and 4\n";

    //STEP 4 : LIST DEDUPLICATION + MAPPING ====================================
    {
        if (build_parameters.verbose > 0) std::cerr << "Step 4: list deduplication and mapping\n";
        start_time = std::chrono::high_resolution_clock::now();
        
        typename ColorClasses::builder cbuild(m_filenames.size(), build_parameters.verbose);
        pthash::compact_vector::builder m_map_builder(hf.num_keys(), ceil(log2(hf.num_keys()))+build_parameters.b); //1bit for check
        // TODO: ceil(log2(hf.num_keys())) depends on the number of unique minmer, should depend on the number of distinct colors instead, but should not bug because nb_distinct_colors <= nb_unique_minmers

        if (build_parameters.check) {
            if (build_parameters.verbose > 0) std::cerr << "map/MPHF of size: " << hf.num_keys() << "\n";
            for (std::size_t i = 0; i < hf.num_keys(); ++i) 
                m_map_builder.set(i, (1 << m_map_builder.width()) - 1); //max value with ceil(log2(hf.num_keys())) bits
        }

        uint32_t cid = 0;
        uint32_t cid_with_parity = 0;

        auto itr = final_result.cbegin();
        while(itr != final_result.cend()) {
            auto current_color = (*itr).first;
            cbuild.add_color_set(current_color.data(), current_color.size()); // only one copy of the color in storage
            while(itr != final_result.cend() and (*itr).first == current_color) {
                auto minimizer = (*itr).second;
                auto mp_idx = hf(minimizer);
                // std::cerr << minimizer << " -> " << mp_idx << "\n";
                cid_with_parity = (cid << build_parameters.b) | ( minimizer & ((1UL << build_parameters.b)-1) );
                m_map_builder.set(mp_idx, cid_with_parity);
                ++itr;
            }
            ++cid;
        }
        cbuild.build(m_ccs);
        m_map_builder.build(m_map);
        //compact_vector m_map where m_map[ hf(minmer) ] = color_id 
    }

    std::cerr << "DEBUG Time for dedup + mapping: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::high_resolution_clock::now() - start_time
                 ).count() 
              << " milliseconds\n";

    std::cerr << "DEBUG end mapping\n";
    std::cerr << "check mem breakdown for thr rest \n";

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