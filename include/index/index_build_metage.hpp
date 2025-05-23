#ifndef INDEX_BUILD_METAGE_HPP
#define INDEX_BUILD_METAGE_HPP

namespace kaminari {
namespace minimizer {

#define CLASS_HEADER template <class ColorClasses, class ColorMapper>
#define METHOD_HEADER index<ColorClasses, ColorMapper>


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
METHOD_HEADER::read_file_task(const std::string& file, color_t doc_id, emem_t& result, const build::options_t& build_parameters) {
    // std::size_t total_kmers = 0;
    // std::size_t total_mmers = 0;
    // std::size_t total_minimizers = 0;
    gzFile fp = nullptr;
    kseq_t* seq = nullptr;
    std::vector<::minimizer::record_t> mms_buffer;
    if ((fp = gzopen(file.c_str(), "r")) == NULL)
        throw std::runtime_error("Unable to open input file " + file);
    seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        std::size_t contig_mmer_count;
        [[maybe_unused]] auto contig_kmer_count = ::minimizer::from_string<hash64>(
            seq->seq.s, 
            seq->seq.l, 
            build_parameters.k, 
            build_parameters.m, 
            build_parameters.seed, 
            build_parameters.canonical,
            contig_mmer_count,
            mms_buffer
        );
        // total_kmers += contig_kmer_count;
        // total_mmers += contig_mmer_count;
        // total_minimizers += mms_buffer.size();

        // duplicates removed during merge
        std::vector<minimizer_t> minimizers;
        for (auto r : mms_buffer) {
            result.push_back(
                std::make_pair(r.hash, doc_id));
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
METHOD_HEADER::build_metage(const build::options_t& build_parameters)
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
    color_t doc_id = 0;
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
    for (size_t i = 0; i < build_parameters.nthreads; ++i) {
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

    std::cerr << "Time for parsing: " 
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - start_time
                    ).count() 
                << " milliseconds\n";
    // all emem created, filled, sorted, dumped to disk

    //STEP 2 : MERGE RESULTS ===================================================
    if (build_parameters.verbose > 0) std::cerr << "Step 2: Merging results of parsing\n";
    start_time = std::chrono::high_resolution_clock::now();
    ankerl::unordered_dense::set<minimizer_t> unique_minmers;

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

    minimizer_t smallest;
    while (itr != end){
        smallest = (*itr).first;
        std::pair<std::vector<color_t>, minimizer_t> color_minmer = {{(*itr).second}, smallest};
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

    std::cerr << "Time for merging: " 
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

        std::cerr << "Time for MPHF build: " 
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - start_time
                    ).count() 
                << " milliseconds\n";

        std::cerr << "size of MPHF : " << hf.num_bits()/8 << " Bytes\n";
    }

    //STEP 4 : LIST DEDUPLICATION + MAPPING ====================================
    {
        if (build_parameters.verbose > 0) std::cerr << "Step 4: list deduplication and mapping\n";
        start_time = std::chrono::high_resolution_clock::now();
        
        typename ColorClasses::builder cbuild(m_filenames.size(), build_parameters.verbose);
        kaminari::compact_vector::builder m_map_builder(hf.num_keys(), ceil(log2(hf.num_keys()))+build_parameters.b); //1bit for check
        // TODO: ceil(log2(hf.num_keys())) depends on the number of unique minmer, should depend on the number of distinct colors instead, but should not bug because nb_distinct_colors <= nb_unique_minmers

        color_t cid = 0;
        color_t cid_with_parity = 0;

        uint64_t minimizer;
        uint64_t mp_idx;


        auto itr = final_result.cbegin();
        while(itr != final_result.cend()) {
            auto current_color = (*itr).first;
            cbuild.add_color_set(current_color.data(), current_color.size()); // only one copy of the color in storage
            while(itr != final_result.cend() and (*itr).first == current_color) {
                minimizer = (*itr).second;
                mp_idx = hf(minimizer);
                // std::cerr << minimizer << " -> " << mp_idx << "\n";
                cid_with_parity = (cid << build_parameters.b) | ( minimizer & ((1UL << build_parameters.b)-1) );
                m_map_builder.set(mp_idx, cid_with_parity);
                ++itr;
            }
            ++cid;
        }

        //shrink color mapper because we store cids so we only need log2(cid) bits
        m_map_builder.shrink(ceil(log2(hf.num_keys())) - ceil(log2(cid)));

        cbuild.build(m_ccs);
        m_map_builder.build(m_map);
        //compact_vector m_map where m_map[ hf(minmer) ] = color_id 
    }

    std::cerr << "Time for dedup + mapping: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::high_resolution_clock::now() - start_time
                 ).count() 
              << " milliseconds\n";

    if (build_parameters.verbose > 0) {
        std::cerr << "Number of ids: " << m_ccs.num_docs() << "\n";
        std::cerr << "Number of colors (lists of ids):" << m_ccs.num_color_classes() << "\n";
    } 
}




#undef CLASS_HEADER
#undef METHOD_HEADER

} // namespace minimizer
} // namespace kaminari

#endif // INDEX_BUILD_METAGE_HPP