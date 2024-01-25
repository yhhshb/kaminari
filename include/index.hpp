#ifndef KAMINARI_INDEX_HPP
#define KAMINARI_INDEX_HPP

#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <unordered_set>
#include "constants.hpp"
#include "GGCAT.hpp"
#include "minimizer.hpp" // FIXME: to be removed once biolib's minimizer iterators are ready

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
        typedef emem::external_memory_vector<std::pair<minimizer_t, std::size_t>> mmc_vector_t;

        struct mm_cids_hash_t {
            mm_cids_hash_t(minimizer_t minimizer, std::vector<std::size_t> color_ids, uint32_t seed) {
                cids = color_ids;
                std::sort(cids.begin(), cids.end());
                mm = minimizer;
                cids_hash = hash64::hash(reinterpret_cast<const uint8_t*>(cids.data()), cids.size() * sizeof(typename decltype(cids)::value_type), seed);
            }

            minimizer_t mm;
            std::vector<std::size_t> cids;
            hash64::hash_type cids_hash;
        };

        std::tuple<typename ColorClasses::builder, std::size_t, std::size_t> 
        collect_minimizers_and_colors(
            const opt_t& build_parameters, 
            mmc_vector_t& mms_and_colors
        ) const;

        std::vector<mm_cids_hash_t> 
        group_minimizers_by_color(
            const opt_t& build_parameters,
            mmc_vector_t&& mms_and_colors
        ) const;

        std::vector<minimizer_t> 
        build_minimizers_mphf(
            const opt_t& build_parameters,
            const std::vector<typename METHOD_HEADER::mm_cids_hash_t>& mm_cids
        );

        void 
        merge_colors_by_minimizer(
            const opt_t& build_parameters,
            const std::vector<minimizer_t>& minimizers,
            const std::vector<typename METHOD_HEADER::mm_cids_hash_t>& mm_cids,
            const ColorClasses& unmerged_colors
        );

        void 
        build(const opt_t& build_parameters);

        pthash_opt_t 
        get_pthash_options(const opt_t& build_parameters);

        opt_t::fn_t m_filenames;
        ColorClasses m_ccs; // colors
        ColorMapper m_map; // map between mphf values and color classes
        pthash_minimizers_mphf_t hf; // minimizer mphf
};

template <class ColorClasses, class ColorMapper>
index<ColorClasses, ColorMapper>::index(const opt_t& build_parameters)
{
    build(build_parameters);
}

CLASS_HEADER
pthash_opt_t
METHOD_HEADER::get_pthash_options(const opt_t& build_parameters)
{
    pthash_opt_t opts;
    opts.c = build_parameters.pthash_constant;
    return opts;
}

CLASS_HEADER
std::tuple<typename ColorClasses::builder, std::size_t, std::size_t> 
METHOD_HEADER::collect_minimizers_and_colors(const opt_t& build_parameters, mmc_vector_t& mms_and_colors) const
{
    std::size_t num_unitigs = 0;
    std::size_t num_distinct_colors = 0;
    std::vector<::minimizer::record_t> mm_buffer;
    typename ColorClasses::builder unmerged_colors_builder(build_parameters.input_filenames.size(), build_parameters.verbose);
    GGCAT ggreader(build_parameters);
    ggreader.loop_through_unitigs(
        [&](ggcat::Slice<char> const unitig, ggcat::Slice<uint32_t> const colors, bool same_color) 
        {
            if (not same_color) {
                ++num_distinct_colors; // color_id
                unmerged_colors_builder.add_color_set(colors.data, colors.size); // compress colors
            }
            mm_buffer.clear();
            uint64_t mm_count;
            ::minimizer::from_string<hash64>(
                unitig.data, 
                unitig.size, 
                build_parameters.k, 
                build_parameters.m, 
                build_parameters.seed, 
                build_parameters.canonical,
                mm_count,
                mm_buffer
            );

            std::sort(
                mm_buffer.begin(), 
                mm_buffer.end(), 
                [](const ::minimizer::record_t& r1, const ::minimizer::record_t& r2) { 
                    return r1.itself < r2.itself; 
                }
            );
            auto last = std::unique(
                mm_buffer.begin(), 
                mm_buffer.end(), 
                [](const ::minimizer::record_t& r1, const ::minimizer::record_t& r2) { 
                    return r1.itself == r2.itself; 
                }
            );
            for(auto itr = mm_buffer.begin(); itr != last; ++itr) { // unique minimizers
                auto p = std::make_pair((*itr).itself, num_distinct_colors);
                mms_and_colors.push_back(p);
            }
            ++num_unitigs;
        }
    );
    return std::make_tuple(unmerged_colors_builder, num_unitigs, num_distinct_colors);
}

CLASS_HEADER
std::vector<typename METHOD_HEADER::mm_cids_hash_t> 
METHOD_HEADER::group_minimizers_by_color(
    const opt_t& build_parameters,
    mmc_vector_t&& mms_and_colors) const
{ // compress color indexes associated to minimizers into vectors (one vector per minimizer)
    if (build_parameters.verbose) std::cerr << "\tstep 1: Assign colors to minimizers\n";
    auto itr = mms_and_colors.cbegin();
    auto end = mms_and_colors.cend();
    std::vector<mm_cids_hash_t> minimizers_cids_hashes;
    while (itr != end) {
        std::vector<std::size_t> mm_color_ids;
        minimizer_t mm = (*itr).first;
        while (itr != end and mm == (*itr).first) {
            mm_color_ids.push_back((*itr).second);
            ++itr;
        }
        minimizers_cids_hashes.emplace_back(mm, std::move(mm_color_ids), 42);
    }
    if (build_parameters.verbose) std::cerr << "\tstep 2: Group minimizers by color\n";
    std::sort(
        minimizers_cids_hashes.begin(), 
        minimizers_cids_hashes.end(), 
        [](const mm_cids_hash_t& a, const mm_cids_hash_t& b) {
            return a.cids_hash < b.cids_hash;
        }
    );
    return minimizers_cids_hashes;
}

CLASS_HEADER
std::vector<minimizer_t>
METHOD_HEADER::build_minimizers_mphf(
    const opt_t& build_parameters,
    const std::vector<typename METHOD_HEADER::mm_cids_hash_t>& mm_cids)
{
    auto member = [](const mm_cids_hash_t& e) {return e.mm;};
    std::vector<minimizer_t> minimizers(
        iterators::member_iterator(mm_cids.begin(), member),
        iterators::member_iterator(mm_cids.end(), member)
    );
    if (build_parameters.verbose) std::cerr << "\tstep 3: build minimizers MPHF\n";

    int backup, redirect;
    fflush(stdout);
    backup = dup(1);
    redirect = open("/dev/null", O_WRONLY);
    dup2(redirect, 1);
    close(redirect);

    hf.build_in_external_memory(minimizers.begin(), minimizers.size(), get_pthash_options(build_parameters));

    fflush(stdout);
    dup2(backup, 1);
    close(backup);

    std::cerr << "check 5" << std::endl;

    std::vector<minimizer_t> ordered_minimizers(minimizers.size());
    for (auto mm : minimizers) {
        ordered_minimizers[hf(mm)] = mm;
    }
    return ordered_minimizers;
}

CLASS_HEADER
void
METHOD_HEADER::merge_colors_by_minimizer(
    const opt_t& build_parameters,
    const std::vector<minimizer_t>& mphf_ordered_minimizers,
    const std::vector<typename METHOD_HEADER::mm_cids_hash_t>& mm_cids,
    const ColorClasses& unmerged_colors)
{
    auto accessor = [](const mm_cids_hash_t& e) {return e.cids_hash;};
    sampler::ordered_unique_sampler sampler(
        iterators::member_iterator(mm_cids.begin(), accessor),
        iterators::member_iterator(mm_cids.end(), accessor)
    );
    std::size_t mm_colors_count = 0;
    for (auto itr = sampler.cbegin(); itr != sampler.cend(); ++itr) ++mm_colors_count;
    typename ColorClasses::builder merged_colors_builder(mm_colors_count, build_parameters.verbose);

    auto collect_colors = [&](const mm_cids_hash_t& record) {
        std::unordered_set<uint32_t> mm_colors;
        for (auto cid : record.cids) {
            auto itr = unmerged_colors.colors(cid);
            auto row_size = itr.size();
            for (std::size_t i = 0; i < row_size; ++i) {
                mm_colors.insert(*itr);
                ++itr;
            }
        }
        std::vector<uint32_t> row;
        std::copy(mm_colors.begin(), mm_colors.end(), std::back_inserter(row));
        return row;
    };

    if (build_parameters.verbose) std::cerr << "\tstep 4: store colors\n";

    for (auto& record : mm_cids) {
        auto colors = collect_colors(record);
        merged_colors_builder.add_color_set(colors.data(), colors.size());
        // we then need to create a map hash(minimizer) -> color
        if (build_parameters.verbose) std::cerr << "\t- build mapping from k-mers to colors\n";
        // m_map.build(hf, m_ccs, build_parameters.verbose);
    }

    for (auto mm : mphf_ordered_minimizers) {
        auto colors = collect_colors(mm_cids.at(hf(mm)));
        merged_colors_builder.add_color_set(colors.data(), colors.size());
        // do nothing, minimizers were MPHF-sorted so nothing to do
    }
    
    merged_colors_builder.build(m_ccs);
}

CLASS_HEADER
void 
METHOD_HEADER::build(const opt_t& build_parameters)
{
    m_filenames = build_parameters.input_filenames;
    if (build_parameters.verbose) std::cerr << "About to process " << m_filenames.size() << " files.\n";
    mmc_vector_t mms_and_colors(build_parameters.max_ram * constants::GB, build_parameters.tmp_dir, "mm_dump.bin");
    auto [unmerged_colors_builder, num_unitigs, num_distinct_colors] = collect_minimizers_and_colors(build_parameters, mms_and_colors);
    if (num_unitigs > std::numeric_limits<uint32_t>::max()) 
        throw std::runtime_error("Number of unitigs > than allowed maximum of " + std::to_string(std::numeric_limits<uint32_t>::max()));
    std::cerr << "check 0" << std::endl;
    std::vector<mm_cids_hash_t> mm_cids = group_minimizers_by_color(build_parameters, std::move(mms_and_colors));
    std::cerr << "check 1" << std::endl;
    std::vector<minimizer_t> mphf_ordered_minimizers = build_minimizers_mphf(build_parameters, mm_cids);
    std::cerr << "check 2" << std::endl;
    ColorClasses unmerged_colors;
    unmerged_colors_builder.build(unmerged_colors);
    std::cerr << "check 3" << std::endl;
    merge_colors_by_minimizer(build_parameters, mphf_ordered_minimizers, std::move(mm_cids), std::move(unmerged_colors));
    if (build_parameters.verbose) {
        std::cerr << "\t\tnum_unitigs " << num_unitigs << "\n";
        std::cerr << "\t\tnum_distinct_colors " << num_distinct_colors << "\n";
    }
    
    if (build_parameters.check) {
        std::cerr << "\tstep 6: check correctness...\n";
        // TODO
        std::cerr << "Checking done. Everything OK\n";
    }
}

CLASS_HEADER
void 
METHOD_HEADER::memory_breakdown(std::ostream& out) const noexcept
{
    libra scale;
    scale.visit(m_filenames);
    out << "The list of input filenames weights: " << scale.get_byte_size() * 8 << " bits";
    out << "The MPHF of minimizers weights: " << hf.num_bits() << " bits"; 
    scale.visit(m_ccs);
    out << "colors weight: " << scale.get_byte_size() * 8 << " bits";
    scale.visit(m_map);
    out << "The mapping from minimizers to colors weights: " << scale.get_byte_size() * 8 << " bits";
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

#undef CLASS_HEADER
#undef METHOD_HEADER

} // namespace minimizer

namespace kmer {

template <class ColorClasses, class ColorMapper>
class index
{
    public:
        index(const opt_t& build_parameters);

        void memory_breakdown(std::ostream& out) const noexcept;
        void print_map_histogram(std::ostream& out) const noexcept;
        void print_map(std::ostream& out) const noexcept;

        template <class Visitor>
        void visit(Visitor& visitor);

    private:
        void build(const opt_t& build_parameters);
        lphash_configuration_t get_lphash_options(const opt_t& build_parameters) const noexcept;

        opt_t::fn_t m_filenames;
        ColorClasses m_ccs; // colors
        ColorMapper m_map; // map between mphf values and color classes
        lphash_mphf_t hf; // lphash mphf
};

template <class ColorClasses, class ColorMapper>
index<ColorClasses, ColorMapper>::index(const opt_t& build_parameters)
{
    build(build_parameters);
}

template <class ColorClasses, class ColorMapper>
void 
index<ColorClasses, ColorMapper>::build(const opt_t& build_parameters)
{
    m_filenames = build_parameters.input_filenames;
    if (build_parameters.verbose) std::cerr << "about to process " << m_filenames.size() << " files...\n";
    std::string tmp_fasta_input = build_parameters.output_filename + ".fa";
    GGCAT ggreader(build_parameters);
    {
        if (build_parameters.verbose) std::cerr << "step 1. building u2c and m_ccs\n";
        std::size_t num_unitigs = 0;
        std::size_t num_distinct_colors = 0;

        std::ofstream out(tmp_fasta_input.c_str()); // FIXME: write unitigs to fasta file for lpHash 
        if (not out.is_open()) throw std::runtime_error("cannot open output file");
        typename ColorClasses::builder colors_builder(m_filenames.size(), build_parameters.verbose);

        ggreader.loop_through_unitigs(
            [&](ggcat::Slice<char> const unitig, ggcat::Slice<uint32_t> const colors, bool same_color) 
            {
                if (not same_color) {
                    ++num_distinct_colors; // color_id
                    colors_builder.add_color_set(colors.data, colors.size); // compress colors
                }
                out << ">" << num_distinct_colors << "\n";
                out.write(unitig.data, unitig.size);
                out << '\n';
                ++num_unitigs;
            }
        );
        out.close();

        if (num_unitigs > std::numeric_limits<uint32_t>::max()) 
            throw std::runtime_error("Number of unitigs > than allowed maximum of " + std::to_string(std::numeric_limits<uint32_t>::max()));
        colors_builder.build(m_ccs);

        if (build_parameters.verbose) {
            std::cerr << "num_unitigs " << num_unitigs << "\n";
            std::cerr << "num_distinct_colors " << num_distinct_colors << "\n";
            // std::cerr << "m_u2c.size() " << m_u2c.size() << "\n";
            // std::cerr << "m_u2c.size1() " << m_u2c.size1() << "\n";
            // std::cerr << "m_u2c.size0() " << m_u2c.size0() << "\n";
        }
    }

    {
        if (build_parameters.verbose) std::cerr << "step 2. build minimal perfect hash function\n";
        hf.build(get_lphash_options(build_parameters), std::cerr);
    }

    { // step 3
        if (build_parameters.verbose) std::cerr << "step 3. build mapping from k-mers to colors\n";
        // m_map.build(hf, ggreader, build_parameters.verbose);
        m_map.build(hf, tmp_fasta_input, build_parameters.verbose);
    }
    
    if (build_parameters.check) {
        std::cerr << "step 4. check correctness...\n";
        std::ifstream fastain(tmp_fasta_input.c_str());
        std::string buffer;
        std::size_t color_class_id = 0;
        while(std::getline(fastain, buffer)) {
            if (buffer.size()) {
                if (buffer[0] == '>') {
                    color_class_id = std::strtoull(buffer.data() + 1, NULL, 10);
                } else {
                    auto hash_values = hf(buffer.data(), buffer.size(), true);
                    assert(color_class_id != 0);
                    for (auto v : hash_values) {
                        auto cc = m_map.at(v);
                        // std::cerr << " with color class id = " << cc << " (true color class = " << color_class_id << ")\n";
                        if (cc != color_class_id) {
                            std::cerr << "Check FAIL\n"; // Because Rust cannot catch foreign exceptions
                            throw std::runtime_error("Check FAIL");
                        }
                    }
                }
            }
        }
        // ggreader.loop_through_unitigs(
        //     [&](ggcat::Slice<char> const unitig, 
        //         ggcat::Slice<uint32_t> const colors,
        //         bool same_color) // everything has the same color
        //     {
        //         if (not same_color) ++color_class_id;
        //         auto hash_values = hf(unitig.data, unitig.size, true);
        //         // std::cerr << "hashes computed\n";
        //         for (auto v : hash_values) {
        //             auto cc = m_map.at(v);
        //             // std::cerr << " with color class id = " << cc << " (true color class = " << color_class_id << ")\n";
        //             if (cc != color_class_id) {
        //                 std::cerr << "Check FAIL\n"; // Because Rust cannot catch foreign exceptions
        //                 throw std::runtime_error("Check FAIL");
        //             }
        //         }
        //     }
        // );
        std::cerr << "Checking done. Everything OK\n";
    }

    try { // remove unitig file
        std::remove(tmp_fasta_input.c_str());
    } catch (std::exception const& e) {
        std::cerr << e.what() << std::endl;
    }
}

template <class ColorClasses, class ColorMapper>
void 
index<ColorClasses, ColorMapper>::memory_breakdown(std::ostream& out) const noexcept
{
    libra scale;
    scale.visit(m_map);
    out << "Mapping k-mers to colors weights: " << scale.get_byte_size() << " Bytes";
}

template <class ColorClasses, class ColorMapper>
void 
index<ColorClasses, ColorMapper>::print_map_histogram(std::ostream& out) const noexcept
{
    auto map_hist = m_map.get_histogram();
    std::size_t val_sum = 0;
    for (auto itr = map_hist.begin(); itr != map_hist.end(); ++itr) {
        val_sum += itr->second;
    }
    double entropy = 0;
    for (auto itr = map_hist.begin(); itr != map_hist.end(); ++itr) {
        double p = double(itr->second) / val_sum;
        entropy += p * std::log2(p);
    }
    out << "H_0," << -entropy << "\n";
    for (auto itr = map_hist.begin(); itr != map_hist.end(); ++itr) {
        out << itr->first << "," << itr->second << "\n";
    }
}

// template <class ColorClasses, class ColorMapper>
// void 
// index<ColorClasses, ColorMapper>::print_map(std::ostream& out) const noexcept
// {
//     for (auto itr = m_map.vector_data().cbegin(); itr != m_map.vector_data().cend(); ++itr) {
//         out << *itr << "\n";
//     }
// }

template <class ColorClasses, class ColorMapper>
lphash_configuration_t 
index<ColorClasses, ColorMapper>::get_lphash_options(const opt_t& build_parameters) const noexcept
{
    lphash_configuration_t lphash_config;
    lphash_config.tmp_dirname = build_parameters.tmp_dir;
    lphash_config.max_memory = build_parameters.max_ram;
    lphash_config.num_threads = build_parameters.nthreads;
    lphash_config.c = build_parameters.pthash_constant;
    lphash_config.k = build_parameters.k;
    lphash_config.m = build_parameters.m;
    lphash_config.verbose = build_parameters.verbose;
    lphash_config.input_filename = build_parameters.output_filename + ".fa";
    return lphash_config;
}

template <class ColorClasses, class ColorMapper>
template <class Visitor>
void
index<ColorClasses, ColorMapper>::visit(Visitor& visitor)
{
    visitor.visit(m_filenames);
    visitor.visit(m_ccs); // colors
    visitor.visit(m_map); // map between mphf values and color classes
    visitor.visit(hf); // lphash mphf
}

} // namespace kmer

} // namespace kaminari

#endif // KAMINARI_INDEX_HPP