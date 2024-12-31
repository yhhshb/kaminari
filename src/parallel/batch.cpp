#include "../../include/index.hpp"
#include "../../include/hybrid.hpp"

namespace kaminari {
namespace minimizer {

#define CLASS_HEADER template <class ColorClasses, class ColorMapper>
#define METHOD_HEADER index<ColorClasses, ColorMapper>

template class index<kaminari::color_classes::hybrid, pthash::compact_vector>;



/*
STEP 1 PROCESS FILES : PARSE FILES, THEN MERGE BATCH (ONLY 1)
=========================================================
*/
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
    {
        std::lock_guard<std::mutex> lock(task_stack_mutex);
        for (uint32_t docid = start; docid < end; ++docid) {
            //create read_file_task for each file of batch
            task_stack.push( //Function(f, desc)
                {
                    [this, &start, allocated_ram, docid, &d1_results_storage, &d1_results_storage_mutex, &running_task, &cv, &build_parameters, &debug_cerr_mutex]() {
                        ++running_task;

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
METHOD_HEADER::run_first_batch(
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
    //2 steps : After init_batch, files parsed and stored into emem vector, first merging these vectors then results (emem vector of pairs minmer/[docids]) are iteratively merged until only one remains
    
    std::deque<depthn_result> tmp_results_storage;
    std::mutex tmp_results_storage_mutex;
    uint16_t depth = 1;
    uint64_t ram_for_1_depth = build_parameters.max_ram * constants::GB / 2;
    uint64_t allocated_ram;
    uint32_t nb_results = d1_results_storage.size();

    if (m_filenames.size() == 1) {
        std::cerr << "TODO: handle case 1 file in index, likely going to crash\n";
    }

    

    //1ST STEP : PROCESS DEPTH 1 RESULTS
    allocated_ram = (nb_results == 1) ? ram_for_1_depth : ram_for_1_depth / (nb_results/2); //divide by the nb of results in next depth
    if (nb_results % 2 == 1) {
        //odd number of files in batch, possibly last batch (if even number of threads)
        if (nb_results == 1) {
            //special case last batch contains only one file, merge it with last batch's result
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

        else {  //odd number (make it even by merging 3 firsts together)
            std::lock_guard<std::mutex> lock(task_stack_mutex);
            task_stack.push(
                {
                    [this, &batch_id, &depth, allocated_ram, &d1_results_storage, &d1_results_storage_mutex, &tmp_results_storage, &tmp_results_storage_mutex, &build_parameters, &cv, &running_task, &debug_cerr_mutex]() {
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
            nb_results = nb_results - 3; 
        }
    } //1ST STEP : end case odd number of files in batch

    { //1ST STEP : even number of files in batch, merge them 2 by 2
        std::lock_guard<std::mutex> lock(task_stack_mutex);
        while (nb_results != 0) {
            task_stack.push(
                {
                    [this, &batch_id, &depth, allocated_ram, &d1_results_storage, &d1_results_storage_mutex, &tmp_results_storage, &tmp_results_storage_mutex, &build_parameters, &cv, &running_task, &debug_cerr_mutex]() {
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
            nb_results = nb_results - 2;
        }
    }

    //1ST STEP : all merging tasks added to stack, notify threads
    cv.notify_all();

    //1ST STEP : wait for all tasks to finish
    {
        std::unique_lock<std::mutex> lock(task_stack_mutex);
        cv.wait(lock, [&]() {
            return running_task.load() == 0 && task_stack.empty();
        });
    }
    
    //2ND STEP : MERGE DEPTH N RESULTS, FINISH FIRST BATCH
    depth++;
    nb_results = tmp_results_storage.size();
    while (nb_results != 1){
        allocated_ram = (nb_results == 1) ? ram_for_1_depth : ram_for_1_depth / (nb_results/2);
        //iteratively go deeper into batch until everything merged 
        {
            std::lock_guard<std::mutex> lock(task_stack_mutex);
            if (nb_results % 2 == 1) {
                task_stack.push(
                    {
                        [this, &batch_id, &depth, allocated_ram, &tmp_results_storage, &tmp_results_storage_mutex, &build_parameters, &cv, &running_task, &debug_cerr_mutex]() {
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
            }//2ND STEP : end case odd number of results in batch

            while (nb_results != 0) {
                ///2ND STEP : even number of results, merge them 2 by 2
                task_stack.push(
                    {
                        [this, &batch_id, &depth, allocated_ram, &tmp_results_storage, &tmp_results_storage_mutex, &build_parameters, &cv, &running_task, &debug_cerr_mutex]() {
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
        } //2ND STEP : end lock on stack
        cv.notify_all(); // 2ND STEP : Notify threads about new tasks

        //2ND STEP : wait for all tasks to finish
        {
            std::unique_lock<std::mutex> lock(task_stack_mutex);
            cv.wait(lock, [&]() {
                return running_task.load() == 0 && task_stack.empty();
            });
        }
            

        nb_results = tmp_results_storage.size();
        depth++;
        //2ND STEP : start again until only 1 result left
    }

    //all results merged, transfer to common results_storage
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
    //Calling init_batch() (parse) and run_first_batch() (merge) for each batch of files
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











/*
STEP 2 PROCESS TREE : MERGE DEPTH N RESULTS UNTIL FINAL RESULT
===========================================================
*/
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
                            [this, &batch_id, &depth, allocated_ram, &tmp_results_storage, &tmp_results_storage_mutex, &build_parameters, &cv, &running_task, &debug_cerr_mutex]() {
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
                        [this, &batch_id, &depth, allocated_ram, &tmp_results_storage, &tmp_results_storage_mutex, &build_parameters, &cv, &running_task, &debug_cerr_mutex]() {
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
                        [this, &batch_id, &depth, allocated_ram, &tmp_results_storage, &tmp_results_storage_mutex, &build_parameters, &cv, &running_task, &debug_cerr_mutex]() {
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



#undef CLASS_HEADER
#undef METHOD_HEADER

} // namespace minimizer
} // namespace kaminari