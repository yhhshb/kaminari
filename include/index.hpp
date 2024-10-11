#ifndef KAMINARI_INDEX_HPP
#define KAMINARI_INDEX_HPP

#include <iostream>
#include "constants.hpp"
#include "utils.hpp"
#include "build_options.hpp"
#include "query_options.hpp"
#include "minimizer.hpp"
#include "../bundled/pthash/include/pthash.hpp"
#include "../bundled/biolib/bundled/prettyprint.hpp"
#include "../bundled/biolib/include/bit_vector.hpp"
#include "../bundled/biolib/include/elias_fano.hpp"

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

        index();
        index(const build::options_t& build_parameters);

        std::vector<color_t> query_full_intersection(char const * const q, const std::size_t l, std::size_t verbosity_level = 0) const noexcept;

        std::vector<color_t> query_full_intersection(const std::string& q, std::size_t verbosity_level = 0) const noexcept {return query_full_intersection(q.c_str(), q.length(), verbosity_level);}


        std::vector<scored_id> query_union_threshold(char const * const q, const std::size_t l, float threshold_ratio, std::size_t verbosity_level = 0) const noexcept;

        std::vector<scored_id> query_union_threshold(const std::string& q, float threshold_ratio, std::size_t verbosity_level = 0) const noexcept {return query_union_threshold(q.c_str(), q.length(), threshold_ratio, verbosity_level);}

        std::vector<scored_id> ranking_query_union_threshold(char const * const q, const std::size_t l, float threshold_ratio, std::size_t verbosity_level = 0) const noexcept;

        std::vector<scored_id> ranking_query_union_threshold(const std::string& q, float threshold_ratio, std::size_t verbosity_level = 0) const noexcept {return ranking_query_union_threshold(q.c_str(), q.length(), threshold_ratio, verbosity_level);}
        
        void memory_breakdown(std::ostream& out) const noexcept;
        
        template <class Visitor>
        void visit(Visitor& visitor);

    private:
        typedef pthash::build_configuration pthash_opt_t;
        typedef pthash::single_phf<pthash::murmurhash2_64, pthash::dictionary_dictionary, true> pthash_minimizers_mphf_t;
        typedef std::pair<minimizer_t, color_t> mpc_t;
        typedef emem::external_memory_vector<mpc_t> mmc_vector_t;

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
        
        std::vector<color_t> full_dense_intersection(std::vector<typename ColorClasses::row_accessor>&& color_id_itrs, std::size_t verbosity_level) const noexcept;
        std::vector<color_t> full_mixed_intersection(std::vector<typename ColorClasses::row_accessor>&& color_id_itrs, std::size_t verbosity_level) const noexcept;

        std::vector<color_t> ranking_dense_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold, uint64_t nb_kmers, std::size_t verbosity_level) const noexcept;
        std::vector<scored_id> ranking_mixed_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold, std::size_t verbosity_level) const noexcept;

        std::vector<color_t> union_dense_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold, uint64_t nb_kmers, std::size_t verbosity_level) const noexcept;
        std::vector<scored_id> union_mixed_intersection(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold, std::size_t verbosity_level) const noexcept;
        
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
    opts.c = build_parameters.pthash_constant;
    opts.seed = build_parameters.seed;
    opts.ram = build_parameters.max_ram * constants::GB;
    opts.num_threads = build_parameters.nthreads;
    opts.tmp_dir = build_parameters.tmp_dir;
    opts.verbose_output = build_parameters.verbose;
    opts.minimal_output = true;
    opts.alpha = 0.94;
    return opts;
}


CLASS_HEADER
void
METHOD_HEADER::build(const build::options_t& build_parameters)
{
    std::size_t run_id = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    mmc_vector_t tmp_sorted_storage(
        build_parameters.max_ram * constants::GB, 
        build_parameters.tmp_dir, 
        utils::get_tmp_filename("", "minimizer_unitig_id", run_id)
    );
    std::unordered_map<uint64_t, std::vector<color_t>> kmer_color_check;
    {
        if (build_parameters.verbose > 0) std::cerr << "Step 1: reading files\n";
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
                if (build_parameters.verbose > 2) {
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
    { // aggregate colors into lists for each minimizer (and build the MPHF while doing so)
        if (build_parameters.verbose > 0) std::cerr << "Step 2: aggregating colors\n";
        emem::external_memory_vector<minimizer_t, false> unique_minimizers(
            build_parameters.max_ram * constants::GB, 
            build_parameters.tmp_dir, 
            utils::get_tmp_filename("", "unique_minimizers", run_id)
        );
        
        auto itr = tmp_sorted_storage.cbegin();
        while (itr != tmp_sorted_storage.cend()) { // build list of colors for each minimizer
            auto current = (*itr).first;
            std::vector<color_t> ids;
            while (itr != tmp_sorted_storage.cend() and (*itr).first == current) {
                ids.push_back((*itr).second);
                ++itr;
            }
            auto last = std::unique(ids.begin(), ids.end());
            ids.erase(last, ids.end());
            sorted_color_lists.push_back(std::make_pair(std::move(ids), current)); // save ([ids], minimizer) to disk
            unique_minimizers.push_back(current);
        }
        if (build_parameters.verbose > 0) std::cerr << "Step 3: building the MPHF for " << unique_minimizers.size() << " minimizers\n";
        auto pt_itr = pthash_input_iterator<decltype(unique_minimizers)::const_iterator>(unique_minimizers.cbegin());
        int backup, redirect;
        fflush(stdout);
        backup = dup(1);
        redirect = open("/dev/null", O_WRONLY);
        dup2(redirect, 1);
        close(redirect);

        hf.build_in_external_memory(pt_itr, unique_minimizers.size(), get_pthash_options(build_parameters));

        fflush(stdout);
        dup2(backup, 1);
        close(backup);
        assert(hf.num_keys() == unique_minimizers.size());
    }

    {
        if (build_parameters.verbose > 0) std::cerr << "Step 4: list deduplication and mapping\n";
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
                    if (check_m_map[mp_idx] != (1 << m_map_builder.width()) - 1) throw std::runtime_error("[check fail] reassigning id of unique minimizer (the minimizer is not unique)");
                }
                m_map_builder.set(mp_idx, cid);
                ++itr;
            }
            ++cid;
        }
        cbuild.build(m_ccs);
        m_map_builder.build(m_map);
    }

    if (build_parameters.verbose > 0) {
        std::cerr << "Number of ids: " << m_ccs.num_docs() << "\n";
        std::cerr << "Number of colors (lists of ids):" << m_ccs.num_color_classes() << "\n";
    }

    /* commented because use query union threshold with scored_id so doesnt compile 
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
                        auto color = query_union_threshold(&seq->seq.s[i], build_parameters.k, 0);
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

#undef CLASS_HEADER
#undef METHOD_HEADER

} // namespace minimizer
} // namespace kaminari

#endif // KAMINARI_INDEX_HPP