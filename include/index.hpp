#ifndef KAMINARI_INDEX_HPP
#define KAMINARI_INDEX_HPP

#include <string>
#include <vector>
#include <iostream>
#include "constants.hpp"
#include "GGCAT.hpp"

namespace kaminari {

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
                // try {
                    if (not same_color) {
                        ++num_distinct_colors; // color_id
                        colors_builder.add_color_set(colors.data, colors.size); // compress colors
                    }
                    out << ">\n";
                    out.write(unitig.data, unitig.size);
                    out << '\n';
                    num_unitigs += 1;
                // } catch (std::exception const& e) {
                //     std::cerr << e.what() << std::endl;
                //     exit(1);
                // }
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
        try { // remove unitig file
            std::remove(tmp_fasta_input.c_str());
        } catch (std::exception const& e) {
            std::cerr << e.what() << std::endl;
        }
    }

    { // step 3
        if (build_parameters.verbose) std::cerr << "step 3. build mapping from k-mers to colors\n";
        m_map.build(hf, ggreader);
    }

    if (build_parameters.check) {
        std::cerr << "step 4. check correctness...\n";
        std::size_t color_class_id = 0;
        ggreader.loop_through_unitigs(
            [&](ggcat::Slice<char> const unitig, 
                ggcat::Slice<uint32_t> const colors,
                bool same_color) // everything has the same color
            {
                // try {
                    if (not same_color) ++color_class_id;
                    auto hash_values = hf(unitig.data, unitig.size, true);
                    for (auto v : hash_values) if (m_map.at(v) != color_class_id) {
                        throw std::runtime_error("Check FAIL");
                    }
                // } catch (std::excpetion const& e) {
                //     std::cerr << e.what() << std::endl;
                //     exit(1);
                // }
            },
            build_parameters.nthreads);
        std::cerr << "Checking done. Everything OK\n";
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

template <class ColorClasses, class ColorMapper>
void 
index<ColorClasses, ColorMapper>::print_map(std::ostream& out) const noexcept
{
    for (auto itr = m_map.vector_data().cbegin(); itr != m_map.vector_data().cend(); ++itr) {
        out << *itr << "\n";
    }
}

template <class ColorClasses, class ColorMapper>
lphash_configuration_t 
index<ColorClasses, ColorMapper>::get_lphash_options(const opt_t& build_parameters) const noexcept
{
    lphash_configuration_t lphash_config;
    lphash_config.tmp_dirname = build_parameters.tmp_dir;
    lphash_config.max_memory = build_parameters.max_ram;
    lphash_config.num_threads = build_parameters.nthreads;
    lphash_config.c = constants::lphash_c;
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

} // namespace kaminari

#endif // KAMINARI_INDEX_HPP