#ifndef KAMINARI_INDEX_HPP
#define KAMINARI_INDEX_HPP

#include <iostream>
#include "constants.hpp"
#include "minimizer.hpp"

#include "../bundled/biolib/bundled/prettyprint.hpp"

namespace kaminari {
namespace minimizer {

#define CLASS_HEADER template <class ColorClasses, class ColorMapper>
#define METHOD_HEADER index<ColorClasses, ColorMapper>

CLASS_HEADER
class index
{
    public:
        index(const opt_t& build_parameters);

        void memory_breakdown(std::ostream& out) const noexcept;
        
        template <class Visitor>
        void visit(Visitor& visitor);

    private:
        typedef std::pair<minimizer_t, color_t> mpc_t;
        typedef emem::external_memory_vector<mpc_t> mmc_vector_t;

        void build(const opt_t& build_parameters);
        pthash_opt_t get_pthash_options(const opt_t& build_parameters);

        opt_t::fn_t m_filenames;
        ColorClasses m_ccs; // colors
        ColorMapper m_map; // map between mphf values and color classes
        pthash_minimizers_mphf_t hf; // minimizer mphf
};

CLASS_HEADER
METHOD_HEADER::index(const opt_t& build_parameters)
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
    visitor.visit(hf); // lphash mphf
    visitor.visit(m_ccs); // colors
    visitor.visit(m_map); // map between mphf values and color classes
}

CLASS_HEADER
void
METHOD_HEADER::build(const opt_t& build_parameters)
{
    std::size_t run_id = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    m_filenames = build_parameters.input_filenames;
    mmc_vector_t tmp_sorted_storage(
        build_parameters.max_ram * constants::GB, 
        build_parameters.tmp_dir, 
        util::get_tmp_filename("", "minimizer_unitig_id", run_id)
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
                uint64_t contig_mmer_count = 0;
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
        util::get_tmp_filename("", "color_minimizer_list", run_id)
    );
    { // aggregate colors into lists for each minimizer (and build the MPHF while doing so)
        if (build_parameters.verbose > 0) std::cerr << "Step 2: aggregating colors\n";
        emem::external_memory_vector<minimizer_t, false> unique_minimizers(
            build_parameters.max_ram * constants::GB, 
            build_parameters.tmp_dir, 
            util::get_tmp_filename("", "unique_minimizers", run_id)
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
        auto pt_itr = non_standard::pthash_input_iterator(unique_minimizers.cbegin());
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
            const std::vector<uint32_t>& current_color = (*itr).first;
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
}

CLASS_HEADER
pthash_opt_t
METHOD_HEADER::get_pthash_options(const opt_t& build_parameters)
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