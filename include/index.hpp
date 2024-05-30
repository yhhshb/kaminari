#ifndef KAMINARI_INDEX_HPP
#define KAMINARI_INDEX_HPP



#include <iostream>
#include "../bundled/pthash/include/pthash.hpp"
#include "utils.hpp"
#include "build_options.hpp"
#include "query_options.hpp"
#include "constants.hpp"
#include "minimizer.hpp"

#include <unordered_map>
#include "../bundled/biolib/bundled/prettyprint.hpp"
#include "../bundled/biolib/include/bit_vector.hpp"
#include "../bundled/biolib/include/elias_fano.hpp"


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
        std::vector<color_t> query_union_threshold(char const * const q, const std::size_t l, float threshold_ratio, std::size_t verbosity_level = 0) const noexcept;
        std::vector<color_t> query_union_threshold(const std::string& q, float threshold_ratio, std::size_t verbosity_level = 0) const noexcept {return query_union_threshold(q.c_str(), q.length(), threshold_ratio);}
        
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
        std::vector<color_t> dense_intersection_victor(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold, uint64_t nb_kmers, std::size_t verbosity_level) const noexcept;
        std::vector<color_t> mixed_intersection_victor(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold, ::size_t verbosity_level) const noexcept;
        
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
                for (auto r : mms_buffer) minimizers.push_back(r.itself);
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
}

CLASS_HEADER
std::vector<typename METHOD_HEADER::color_t> 
METHOD_HEADER::dense_intersection_victor(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold, uint64_t nb_kmers, std::size_t verbosity_level) const noexcept
{
    if (threshold == 0) threshold = 1;
    std::size_t filenames_size = m_filenames.size();
    std::vector<color_t> tmp;

    { // step 1: take the union of complementary sets
        for (auto& itr : color_id_itrs) itr.first.reinit_for_complemented_set_iteration();
        color_t candidate = (
            *std::min_element(
                color_id_itrs.begin(), 
                color_id_itrs.end(),
                [](auto const& x, auto const& y) {
                    return x.first.comp_value() < y.first.comp_value();
                }
            )
        ).first.comp_value();

        uint32_t candidate_count;
        tmp.reserve(filenames_size);
        while (candidate < filenames_size) {
            candidate_count = nb_kmers;
            color_t next_candidate = filenames_size;
            for (uint64_t i = 0; i != color_id_itrs.size(); ++i) {
                if (color_id_itrs.at(i).first.comp_value() == candidate){
                    candidate_count -= color_id_itrs.at(i).second;
                    color_id_itrs[i].first.comp_next();
                } 
                // compute next minimum 
                if (color_id_itrs.at(i).first.comp_value() < next_candidate) {
                    next_candidate = color_id_itrs.at(i).first.comp_value();
                }
            }
            if (candidate_count < threshold) tmp.push_back(candidate);
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
METHOD_HEADER::mixed_intersection_victor(std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>>&& color_id_itrs, uint64_t threshold, std::size_t verbosity_level) const noexcept
{
    if (threshold == 0) threshold = 1; //super low values of opts.threshold_ratio

    std::sort(
        color_id_itrs.begin(),
        color_id_itrs.end(),
        [](auto const& x, auto const& y) { return x.first.size() < y.first.size(); }
    );

    bit::vector<uint64_t> tested(m_filenames.size());
    std::vector<color_t> colors;
    std::size_t vec_size = color_id_itrs.size();
    std::size_t filenames_size = m_filenames.size();
    std::size_t idx;
    std::size_t count;
    bit::ef::array remaining_counts;

    {
        //used to know if worth looking in next colors
        std::vector<uint64_t> ef_build(vec_size+1, 0);
        uint64_t sum = 0;
        for (int i = 1; i <= vec_size; i++) {
            sum += color_id_itrs[vec_size-i].second;
            ef_build[i] += sum;
        }
        ef_build[vec_size] += 1; //acts like a boundary
        remaining_counts = bit::ef::array(ef_build.begin(), vec_size+1, ef_build.back());
    }

    {
        color_t candidate;
        //stop loop when not enough kmers represented in following row_acc
        //to reach threshold 
        for (size_t i = 0; i <= vec_size-1 - remaining_counts.lt_find(threshold); i++){
            if (i>0) color_id_itrs[i].first.reset();
            candidate = color_id_itrs[i].first.value();

            //stop loop when regularly reaching last element
            while (candidate != filenames_size){

                // check if candidate has not been done in other previous row_acc
                //+ need to check if not filenames_size because of tested out_of_range 
                while (candidate != filenames_size && tested.at(candidate)) {
                    color_id_itrs[i].first.next();
                    candidate = color_id_itrs.at(i).first.value();
                }
                // leaving loop when heuristically reaching last element
                if (candidate == filenames_size) {
                    break;
                    //checked every value in this row accessor before next i
                }
                     
                tested.set(candidate); //we chose a never seen candidate, dont check him again

                idx = i+1;
                count = color_id_itrs[i].second;
                
                for (; (idx < vec_size) && (count < threshold); ++idx) {
                    if (i>0) color_id_itrs[idx].first.reset(); //candidate in 2nd row_acc might be lower than candidates in 1st row_acc
                    color_id_itrs[idx].first.next_geq(candidate);
                    color_t val = color_id_itrs.at(idx).first.value();
                    if (val == candidate) {
                        count += color_id_itrs.at(idx).second;
                    } 
                }
                
                if (count >= threshold) {
                    colors.push_back(candidate);
                }

                color_id_itrs.at(i).first.next();
                candidate = color_id_itrs[i].first.value();
            }
        }
    }
    
    return colors;

}

CLASS_HEADER
std::vector<typename METHOD_HEADER::color_t> 
METHOD_HEADER::query_union_threshold(char const * const q, const std::size_t l, float threshold_ratio, std::size_t verbosity_level) const noexcept
{
    if (verbosity_level > 2) std::cerr << "step 1: collect color class ids\n";
    if (verbosity_level > 3) std::cerr << "s : " << q << " l : " << l << "\n"; 
    std::vector<std::pair<std::size_t, uint32_t>> ccids_counts;
    uint64_t contig_kmer_count; 
    { // collect color class ids
        std::size_t contig_mmer_count;
        std::vector<::minimizer::record_t> mms_buffer;
        contig_kmer_count = ::minimizer::from_string<hash64>(q, l, k, m, seed, canonical, contig_mmer_count, mms_buffer);
        for (const auto& record : mms_buffer) { 
            ccids_counts.push_back(std::make_pair(m_map[hf(record.itself)], record.size));
        }
        if (verbosity_level > 3) std::cerr << "query contains " << contig_kmer_count << " k-mers and " << contig_mmer_count << " m-mers\n";
    }
    if (verbosity_level > 3) std::cerr << ccids_counts << "\n";
    if (verbosity_level > 2) std::cerr << "step 2: ids to colors\n";
    std::vector<std::pair<typename ColorClasses::row_accessor, uint32_t>> color_itrs;
    bool all_very_dense = true;


    {
        //sort and unique+erase to remove duplicates, side effect : accumulate counts (n°kmers) for every duplicate
        std::sort(
            ccids_counts.begin(),
            ccids_counts.end(),
            [](auto const& x, auto const& y) { return x.first < y.first; }
        );
        auto last = ccids_counts.erase( kaminari::utils::unique_accumulate( 
                                            ccids_counts.begin(), 
                                            ccids_counts.end() ), 
                                        ccids_counts.end() );

        for (auto itr = ccids_counts.begin(); itr != last; ++itr) { 
            color_itrs.push_back(std::make_pair(m_ccs.colors_at((*itr).first), (*itr).second));
            if (color_itrs.back().first.type() != ColorClasses::row_accessor::complementary_delta_gaps) {
                all_very_dense = false;
            }
        }
    }
    
    if (verbosity_level > 2) {
        std::cerr << "step 3: computing intersections\n";
        if (all_very_dense) std::cerr << "\tcompute dense intersection\n";
        else std::cerr << "\tcompute mixed intersection (for a mix of dense and sparse vectors)\n"; 
    }
    if (color_itrs.empty()) return {};
    if (all_very_dense) return dense_intersection_victor(std::move(color_itrs), contig_kmer_count*threshold_ratio, contig_kmer_count, verbosity_level); // intersect of dense rows
    else return mixed_intersection_victor(std::move(color_itrs), contig_kmer_count*threshold_ratio, verbosity_level); // intersect dense and sparse rows
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
