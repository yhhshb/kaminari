#ifndef KAMINARI_INDEX_HPP
#define KAMINARI_INDEX_HPP

#include <iostream>
#include "../bundled/pthash/include/pthash.hpp"
#include "utils.hpp"
#include "build_options.hpp"
#include "query_options.hpp"
#include "constants.hpp"
#include "minimizer.hpp"

// #include "../bundled/biolib/bundled/prettyprint.hpp"

namespace kaminari {
namespace minimizer {

#define CLASS_HEADER template <class ColorClasses, class ColorMapper>
#define METHOD_HEADER index<ColorClasses, ColorMapper>

CLASS_HEADER
class index
{
    public:
        typedef uint64_t minimizer_t;
        typedef typename ColorClasses::color_t color_t;

        index();
        index(const build::options_t& build_parameters);
        std::vector<color_t> query_full_intersection(char const * const q, const std::size_t l) const noexcept;
        std::vector<color_t> query_full_intersection(const std::string& q) const noexcept {return query_full_intersection(q.c_str(), q.length());}
        std::vector<color_t> query_threshold_union(char const * const q, const std::size_t l, const double threhsold) const noexcept;
        std::vector<color_t> query_threshold_union(const std::string& q, const double threshold) const noexcept {return query_threshold_union(q.c_str(), q.length(), threshold);}
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
        std::vector<color_t> dense_intersection(std::vector<typename ColorClasses::row_accessor>&& color_id_itrs) const noexcept;
        std::vector<color_t> mixed_intersection(std::vector<typename ColorClasses::row_accessor>&& color_id_itrs) const noexcept;
        
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
    scale.visit(m_map);
    out << "The mapping from minimizers to colors weights: " << scale.get_byte_size() << " Bytes\n";
    scale.reset();
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
void
METHOD_HEADER::build(const build::options_t& build_parameters)
{
    std::size_t run_id = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    mmc_vector_t tmp_sorted_storage(
        build_parameters.max_ram * constants::GB, 
        build_parameters.tmp_dir, 
        utils::get_tmp_filename("", "minimizer_unitig_id", run_id)
    );
    {
        if (build_parameters.verbose > 0) std::cerr << "Step 1: reading files\n";
        float fraction = 0.1;
        std::size_t id = 0;
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
                for (auto r : mms_buffer) minimizers.push_back(r.itself);
                std::sort(minimizers.begin(), minimizers.end());
                auto last = std::unique(minimizers.begin(), minimizers.end());
                for (auto itr = minimizers.begin(); itr != last; ++itr) {
                    tmp_sorted_storage.push_back(std::make_pair(*itr, id));
                }
                mms_buffer.clear();
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
    }

    {
        if (build_parameters.verbose > 0) std::cerr << "Step 4: list deduplication and mapping\n";
        typename ColorClasses::builder cbuild(m_filenames.size(), build_parameters.verbose);
        m_map.resize(hf.num_keys());
        uint32_t cid = 0;
        auto itr = sorted_color_lists.cbegin();
        while(itr != sorted_color_lists.cend()) {
            auto current_color = (*itr).first;
            cbuild.add_color_set(current_color.data(), current_color.size()); // only one copy of the color in storage
            while(itr != sorted_color_lists.cend() and (*itr).first == current_color) {
                auto minimizer = (*itr).second;
                m_map[hf(minimizer)] = cid;
                ++itr;
            }
            ++cid;
        }
        cbuild.build(m_ccs);
    }

    if (build_parameters.verbose > 0) {
        std::cerr << "Number of ids: " << m_ccs.num_docs() << "\n";
        std::cerr << "Number of colors (lists of ids):" << m_ccs.num_color_classes() << "\n";
    }
}

CLASS_HEADER
std::vector<typename METHOD_HEADER::color_t> 
METHOD_HEADER::dense_intersection(std::vector<typename ColorClasses::row_accessor>&& color_id_itrs) const noexcept
{
    std::vector<color_t> tmp;
    { // step 1: take the union of complementary sets
        for (auto& itr : color_id_itrs) itr.reinit_for_complemented_set_iteration();
        color_t candidate = (
            *std::min_element(
                color_id_itrs.begin(), 
                color_id_itrs.end(),
                [](auto const& x, auto const& y) {
                    return x.comp_value() < y.comp_value();
                }
            )
        ).comp_value();
        // const uint32_t num_docs = iterators[0].num_docs();

        tmp.reserve(m_filenames.size());
        while (candidate < m_filenames.size()) {
            color_t next_candidate = m_filenames.size();
            for (uint64_t i = 0; i != color_id_itrs.size(); ++i) {
                if (color_id_itrs.at(i).comp_value() == candidate) color_id_itrs[i].comp_next();
                /* compute next minimum */
                if (color_id_itrs.at(i).comp_value() < next_candidate) {
                    next_candidate = color_id_itrs.at(i).comp_value();
                }
            }
            tmp.push_back(candidate);
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
METHOD_HEADER::mixed_intersection(std::vector<typename ColorClasses::row_accessor>&& color_id_itrs) const noexcept
{
    std::sort(
        color_id_itrs.begin(),
        color_id_itrs.end(),
        [](auto const& x, auto const& y) { return x.size() < y.size(); }
    );

    std::vector<color_t> colors;
    {
        // const uint32_t num_docs = iterators[0].num_docs();
        color_t candidate = color_id_itrs.front().value();
        std::size_t i = 1;
        while (candidate < m_filenames.size()) {
            for (; i != color_id_itrs.size(); ++i) {
                color_id_itrs[i].next_geq(candidate);
                color_t val = color_id_itrs.at(i).value();
                if (val != candidate) {
                    candidate = val;
                    i = 0;
                    break;
                }
            }
            if (i == color_id_itrs.size()) {
                colors.push_back(candidate);
                color_id_itrs[0].next();
                candidate = color_id_itrs.at(0).value();
                i = 1;
            }
        }
    }
    return colors;
}

CLASS_HEADER
std::vector<typename METHOD_HEADER::color_t> 
METHOD_HEADER::query_full_intersection(char const * const q, const std::size_t l) const noexcept
{
    std::vector<std::size_t> ccids;
    { // collect color class ids
        std::size_t contig_mmer_count;
        std::vector<::minimizer::record_t> mms_buffer;
        [[maybe_unused]] auto contig_kmer_count = ::minimizer::from_string<hash64>(q, l, k, m, seed, canonical, contig_mmer_count, mms_buffer);
        for (const auto& record : mms_buffer) { 
            ccids.push_back(m_map.at(hf(record.itself)));
        }
    }

    std::vector<typename ColorClasses::row_accessor> color_itrs;
    bool all_very_dense = true;
    {
        auto last = std::unique(ccids.begin(), ccids.end()); // deduplicate color class ids
        for (auto itr = ccids.begin(); itr != last; ++itr) { 
            color_itrs.push_back(m_ccs.colors_at(*itr));

            // IMPROVEMENT 
            // Divide rows based on dense/sparse 
            // 1) apply dense_intersection to the dense dataset 
            // 2) and the complementary iterators from the sparse set <-- Not sure about this
            // make the complementary of the result of (2) and make another intersection for the final result
            if (color_itrs.back().type() != ColorClasses::row_accessor::complementary_delta_gaps) {
                all_very_dense = false;
            }
        }
    }

    if (color_itrs.empty()) return {};
    if (all_very_dense) return dense_intersection(std::move(color_itrs)); // intersect of dense rows
    else return mixed_intersection(std::move(color_itrs)); // intersect dense and sparse rows
}

CLASS_HEADER
std::vector<typename METHOD_HEADER::color_t>
METHOD_HEADER::query_threshold_union(char const * const q, std::size_t l, const double threshold) const noexcept
{
    std::vector<color_t> colors;
    
    return colors;
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

#undef CLASS_HEADER
#undef METHOD_HEADER

} // namespace minimizer
} // namespace kaminari

#endif // KAMINARI_INDEX_HPP