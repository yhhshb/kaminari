#ifndef KAMINARI_INDEX_HPP
#define KAMINARI_INDEX_HPP

#include <iostream>
#include <mutex>
#include <shared_mutex>
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
        void build3(const build::options_t& build_parameters);
        static void process_files(const std::vector<std::string>& files, size_t thread_id, const build::options_t& build_parameters, ankerl::unordered_dense::map<uint64_t, std::vector<color_t>>& local_minmer_to_colors, ankerl::unordered_dense::set<uint64_t>& local_unique_minmers);
        static void process_files2(const std::vector<std::string>& files, size_t thread_id,ankerl::unordered_dense::map<uint64_t, std::vector<color_t>>& minmer_to_colors, std::unordered_map<uint64_t, std::mutex>& locks, const build::options_t& build_parameters);
        void parallel_process_files(ankerl::unordered_dense::map<uint64_t, std::vector<color_t>>& minmer_to_colors, ankerl::unordered_dense::set<uint64_t>& unique_minmers, const build::options_t& build_parameters);
        void parallel_process_files2(ankerl::unordered_dense::map<uint64_t, std::vector<color_t>>& minmer_to_colors, std::unordered_map<uint64_t, std::mutex>& locks, const build::options_t& build_parameters);
        
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
    build3(build_parameters);
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

    /*
    if (build_parameters.check) {
        color_t id = 0;
        gzFile fp = nullptr;
        kseq_t* seq = nullptr;
        std::vector<::minimizer::record_t> mms_buffer;
        for (auto filename : m_filenames) {
            if ((fp = gzopen(filename.c_str(), "r")) == NULL)
                throw std::runtime_error("Unable to open input file " + filename);
            seq = kseq_init(fp);
            while (kseq_read(seq) >= 0) {
                for (std::size_t i = 0; seq->seq.l >= k and i < seq->seq.l - build_parameters.k + 1; ++i) {
                    bool valid_kmer = true;
                    for (auto j = 0; j < k; ++j) {
                        if (::constants::seq_nt4_table[static_cast<uint8_t>(*(seq->seq.s + i + j))] >= 4) valid_kmer = false; 
                    }
                    if (valid_kmer) {
                        auto color = query_union_threshold(&seq->seq.s[i], &seq->seq.s[i].size(), build_parameters);
                        assert(color.size());
                        std::size_t dummy;
                        [[maybe_unused]] auto kmc = ::minimizer::from_string<hash64>(
                            &seq->seq.s[i], 
                            build_parameters.k, 
                            build_parameters.k, 
                            build_parameters.k, 
                            0, 
                            false, 
                            dummy,
                            mms_buffer
                        );
                        assert(mms_buffer.size() <= 1); // there might be N's in the input so 0 valid k-mers
                        if (mms_buffer.size() == 1) {
                            auto color_check = kmer_color_check[mms_buffer.at(0).itself];
                            assert(color_check.size() <= color.size());
                            std::vector<color_t> diff;
                            std::set_symmetric_difference(
                                color.begin(), color.end(),
                                color_check.begin(), color_check.end(),
                                std::back_inserter(diff)
                            );
                            if (diff.size() != (color.size() - color_check.size())) {
                                std::cerr << color << "\n";
                                std::cerr << color_check << "\n";
                                throw std::runtime_error("[check fail] minimizer colors are not a strict superset of the kmer color");
                            }
                        }
                    }
                    mms_buffer.clear();
                }
            }
            if (seq) kseq_destroy(seq);
            gzclose(fp);
            fp = nullptr;
            seq = nullptr;
            ++id;
        }
        std::cerr << "[check PASS] Everything is ok\n";
    }
    */
    
}

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


/* 
================================================================================
BUILD 2 : 1 map for each thread then merge (sequentially for now)
================================================================================
*/

CLASS_HEADER
void
METHOD_HEADER::process_files(const std::vector<std::string>& files,  size_t thread_id, const build::options_t& build_parameters, ankerl::unordered_dense::map<uint64_t, std::vector<color_t>>& local_minmer_to_colors, ankerl::unordered_dense::set<uint64_t>& local_unique_minmers) {
    std::cerr << "\tDEBUG hi thread n°" << thread_id << " here :) \n";

    std::size_t total_kmers = 0;
    std::size_t total_mmers = 0;
    std::size_t total_minimizers = 0;

    gzFile fp = nullptr;
    kseq_t* seq = nullptr;
    std::vector<::minimizer::record_t> mms_buffer;
    color_t id = thread_id * files.size(); // Unique ID per thread chunk

    for (const auto& filename : files) {
        if ((fp = gzopen(filename.c_str(), "r")) == nullptr) {
            throw std::runtime_error("Unable to open input file " + filename);
        }
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
                //avoid duplicates inside a file
                if (local_minmer_to_colors[*itr].empty() || local_minmer_to_colors[*itr].back() != id){
                    local_minmer_to_colors[*itr].push_back(id);
                }
                local_unique_minmers.insert(*itr);
            }
            mms_buffer.clear();
        }
        if (seq) kseq_destroy(seq);
        gzclose(fp);
        fp = nullptr;
        seq = nullptr;
        ++id;
    }
}


CLASS_HEADER
void
METHOD_HEADER::parallel_process_files(ankerl::unordered_dense::map<uint64_t, std::vector<color_t>>& minmer_to_colors, ankerl::unordered_dense::set<uint64_t>& unique_minmers, const build::options_t& build_parameters) {
    auto start_time = std::chrono::high_resolution_clock::now();

    size_t num_threads = build_parameters.nthreads;
    if (num_threads > m_filenames.size()) {
        std::cerr << num_threads << " threads given but only " << m_filenames.size() << " files to process. Using " << m_filenames.size() << " threads for step 1 instead.\n";
        num_threads = m_filenames.size();
    }
    std::vector<std::thread> workers;
    size_t chunk_size = (m_filenames.size() + num_threads - 1) / num_threads; // Round up
    
    std::vector<ankerl::unordered_dense::map<uint64_t, std::vector<color_t>>> thread_minmer_to_colors(num_threads);
    std::vector<ankerl::unordered_dense::set<uint64_t>> thread_unique_minmers(num_threads);

    uint32_t thread_id = 0; 
    auto start = m_filenames.begin();
    while (start + chunk_size < m_filenames.end()){
        std::vector<std::string> chunk_files(start, start + chunk_size);
        workers.push_back(
            std::thread(
                [chunk_files, thread_id, &build_parameters, &thread_minmer_to_colors, &thread_unique_minmers]() {
                    process_files(chunk_files, thread_id, build_parameters, thread_minmer_to_colors[thread_id], thread_unique_minmers[thread_id]);
                }
            )
        );
        start = start + chunk_size;
        thread_id++;
    }
    //last chunk of file, might be empty
    std::vector<std::string> chunk_files(start, m_filenames.end());
    workers.push_back(
        std::thread(
            [chunk_files, thread_id, &build_parameters, &thread_minmer_to_colors, &thread_unique_minmers]() {
                process_files(chunk_files, thread_id, build_parameters, thread_minmer_to_colors[thread_id], thread_unique_minmers[thread_id]);
            }
        )
    );

    
    // Wait for all threads to finish
    for (auto& t : workers) {
        t.join();
    }

    std::cout << "DEBUG Time for all threads to parse: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::high_resolution_clock::now() - start_time
                 ).count() 
              << " milliseconds\n";

    start_time = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < num_threads; ++i) {
        // Merge thread_minmer_to_colors into global_minmer_to_colors
        for (const auto& [key, value] : thread_minmer_to_colors[i]) {
            auto& global_vec = minmer_to_colors[key];
            global_vec.insert(global_vec.end(), value.begin(), value.end());
            std::sort(global_vec.begin(), global_vec.end());
            global_vec.erase(std::unique(global_vec.begin(), global_vec.end()), global_vec.end());
        }

        // Merge thread_unique_minmers into global_unique_minmers
        unique_minmers.insert(thread_unique_minmers[i].begin(), thread_unique_minmers[i].end());
    }
    std::cout << "DEBUG Time to merge threads results " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::high_resolution_clock::now() - start_time
                 ).count() 
              << " milliseconds\n";
}

CLASS_HEADER
void
METHOD_HEADER::build2(const build::options_t& build_parameters)
{
    auto start_time = std::chrono::high_resolution_clock::now();

    ankerl::unordered_dense::map<uint64_t, std::vector<color_t>> minmer_to_colors;
    ankerl::unordered_dense::set<uint64_t> unique_minmers;
    
    //std::size_t run_id = std::chrono::high_resolution_clock::now().time_since_epoch().count();

    parallel_process_files(minmer_to_colors, unique_minmers, build_parameters);

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


/* 
================================================================================
BUILD 3 : 1 shared map and set for all threads
================================================================================
*/

CLASS_HEADER
void
METHOD_HEADER::parallel_process_files2(
    ankerl::unordered_dense::map<uint64_t, std::vector<color_t>>& minmer_to_colors,
    std::unordered_map<uint64_t, std::mutex>& locks,
    const build::options_t& build_parameters
) {
    auto start_time = std::chrono::high_resolution_clock::now();

    size_t num_threads = build_parameters.nthreads;
    if (num_threads > m_filenames.size()) {
        std::cerr << num_threads << " threads given but only " << m_filenames.size() << " files to process. Using " << m_filenames.size() << " threads for step 1 instead.\n";
        num_threads = m_filenames.size();
    }
    std::vector<std::thread> workers;
    size_t chunk_size = (m_filenames.size() + num_threads - 1) / num_threads; // Round up

    uint32_t thread_id = 0; 
    auto start = m_filenames.begin();
    while (start + chunk_size < m_filenames.end()) {
        std::vector<std::string> chunk_files(start, start + chunk_size);
        std::cerr << "DEBUG chunk_files size: " << chunk_files.size() << "\n";
        workers.push_back(
            std::thread(
                [chunk_files, thread_id, &minmer_to_colors, &locks, &build_parameters]() {
                    process_files2(chunk_files, thread_id, minmer_to_colors, locks, build_parameters);
                }
            )
        );
        start = start + chunk_size;
        thread_id++;
    }
    // Process the last chunk
    std::vector<std::string> chunk_files(start, m_filenames.end());
    std::cerr << "DEBUG chunk_files size: " << chunk_files.size() << "\n";
    workers.push_back(
        std::thread(
            [chunk_files, thread_id, &minmer_to_colors, &locks, &build_parameters]() {
                process_files2(chunk_files, thread_id, minmer_to_colors, locks, build_parameters);
            }
        )
    );

    // Wait for all threads to finish
    for (auto& t : workers) {
        t.join();
    }

    std::cout << "DEBUG Time for all threads to parse: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::high_resolution_clock::now() - start_time
                 ).count() 
              << " milliseconds\n";
}

CLASS_HEADER
void
METHOD_HEADER::process_files2(
    const std::vector<std::string>& files,  
    size_t thread_id,
    ankerl::unordered_dense::map<uint64_t, std::vector<color_t>>& minmer_to_colors,
    std::unordered_map<uint64_t, std::mutex>& locks,
    const build::options_t& build_parameters
) {
    std::cerr << "\tDEBUG hi thread n°" << thread_id << " here :) \n";

    std::size_t total_kmers = 0;
    std::size_t total_mmers = 0;
    std::size_t total_minimizers = 0;

    gzFile fp = nullptr;
    kseq_t* seq = nullptr;
    std::vector<::minimizer::record_t> mms_buffer;
    color_t id = thread_id * files.size(); // Unique ID per thread chunk

    for (const auto& filename : files) {
        if ((fp = gzopen(filename.c_str(), "r")) == nullptr) {
            throw std::runtime_error("Unable to open input file " + filename);
        }
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

            std::vector<minimizer_t> minimizers;
            for (auto r : mms_buffer) {
                minimizers.push_back(r.itself);
            }
            std::sort(minimizers.begin(), minimizers.end());
            auto last = std::unique(minimizers.begin(), minimizers.end());
            for (auto itr = minimizers.begin(); itr != last; ++itr) {
                uint64_t minmer = *itr;
                //std::cerr << "\tDEBUG thread" << thread_id << " processing minmer " << minmer << "\n";

                // Ensure a lock exists for this minmer
                {
                    std::lock_guard<std::mutex> lock(locks[minmer]);
                    //std::cerr << "\tDEBUG thread" << thread_id << " locked minmer " << minmer << "\n";
                    if(std::find(minmer_to_colors[minmer].begin(), minmer_to_colors[minmer].end(), id) == minmer_to_colors[minmer].end()) {
                        minmer_to_colors[minmer].push_back(id);
                    } 
                }
            }
            mms_buffer.clear();
        }
        if (seq) kseq_destroy(seq);
        gzclose(fp);
        fp = nullptr;
        seq = nullptr;
        id++;
    }
}

CLASS_HEADER
void
METHOD_HEADER::build3(const build::options_t& build_parameters)
{
    auto start_time = std::chrono::high_resolution_clock::now();

    // Initialize shared data structures
    ankerl::unordered_dense::map<uint64_t, std::vector<color_t>> minmer_to_colors;
    ankerl::unordered_dense::set<uint64_t> unique_minmers;
    std::unordered_map<uint64_t, std::mutex> locks;

    ankerl::unordered_dense::map<uint64_t, uint64_t> first_pass;

    std::size_t total_kmers = 0;
    std::size_t total_mmers = 0;
    std::size_t total_minimizers = 0;

    gzFile fp = nullptr;
    kseq_t* seq = nullptr;
    std::vector<::minimizer::record_t> mms_buffer;

    for (const auto& filename : m_filenames) {
        if ((fp = gzopen(filename.c_str(), "r")) == nullptr) {
            throw std::runtime_error("Unable to open input file " + filename);
        }
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

            std::vector<minimizer_t> minimizers;
            for (auto r : mms_buffer) {
                minimizers.push_back(r.itself);
            }
            std::sort(minimizers.begin(), minimizers.end());
            auto last = std::unique(minimizers.begin(), minimizers.end());
            for (auto itr = minimizers.begin(); itr != last; ++itr) {
                uint64_t minmer = *itr;
                first_pass[minmer]++;
                locks[minmer];
            }
            mms_buffer.clear();
        }
        if (seq) kseq_destroy(seq);
        gzclose(fp);
        fp = nullptr;
        seq = nullptr;
    }

    for (const auto& [minmer, count] : first_pass) {
        minmer_to_colors[minmer].reserve(count);
    }

    std::cout << "DEBUG time first pass: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::high_resolution_clock::now() - start_time
                 ).count() 
              << " milliseconds\n";


    start_time = std::chrono::high_resolution_clock::now();

    parallel_process_files2(minmer_to_colors, locks, build_parameters);

    for (auto& [minmer, colors] : minmer_to_colors) {
        std::sort(colors.begin(), colors.end());
        unique_minmers.insert(minmer);
    }

    std::cout << "DEBUG time second pass: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::high_resolution_clock::now() - start_time
                 ).count() 
              << " milliseconds\n";


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


/*
BUG with buildv2 and t=3, with small ecoli fof, wrong nb of colors


./kaminari_dev build -i fof_ecoli_small.txt -o ecoli_1.kaminari -k 31 -m 19 -t 2 -d ./tmp -g 12 -a -v 1
Step 1: reading files
	total k-mers: 101993239
	total m-mers: 102019320
	total minimizers (with repetitions): 14571168
Step 2: aggregating colors
Step 3: building the MPHF for 2302166 minimizers
Step 4: list deduplication and mapping
m_num_docs: 20
m_sparse_set_threshold_size 5
m_very_dense_set_threshold_size 15
processed 17941 lists
	m_num_total_integers 176827
	total bits for ints = 483992
	total bits per offsets = 133952
	total bits = 617944n	offsets: 0.757531 bits/int
	lists: 2.73709 bits/int
Number of ids: 20
Number of colors (lists of ids):17941
The list of input filenames weights: 1643 Bytes
The MPHF of minimizers weights: 766025 Bytes
Colors weight: 77276 Bytes
The mapping from minimizers to colors weights: 6331000 Bytes

Written 7175978 Bytes*/