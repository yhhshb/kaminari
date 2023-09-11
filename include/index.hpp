#ifndef KAMINARI_INDEX_HPP
#define KAMINARI_INDEX_HPP

#include <string>
#include <vector>
#include "constants.hpp"
#include "GGCAT.hpp"

namespace kaminari{

template <class ColorClasses>
class index
{
    public:
        index(const opt_t& build_parameters);

    private:
        void build(const opt_t& build_parameters);
        opt_t::fn_t m_filenames;
        ranked_bit_vector m_u2c;
        ColorClasses m_ccs;
        // lphash mphf
        // bit vector or not
        // colors
};

template <class ColorClasses>
index<ColorClasses>::index(const opt_t& build_parameters)
{
    build(build_parameters);
}

template <class ColorClasses>
void 
index<ColorClasses>::build(const opt_t& build_parameters)
{
    m_filenames = build_parameters.input_filenames;
    GGCAT ggreader(build_parameters);
    {
        if (build_parameters.verbose) std::cerr << "step 1. building u2c and m_ccs\n";
        std:size_t num_unitigs = 0;
        std::size_t num_distinct_colors = 0;
        bit_vector bvb;  // for m_u2c
        std::ofstream out((build_parameters.output_filename + ".fa").c_str()); // FIXME: write unitigs to fasta file for SSHash 
        if (!out.is_open()) throw std::runtime_error("cannot open output file");

        typename ColorClasses::builder colors_builder(m_filenames.size());

        builder.ggcat->loop_through_unitigs(
            [&](ggcat::Slice<char> const unitig, ggcat::Slice<uint32_t> const colors, bool same_color) 
            {
                try {
                    if (!same_color) {
                        ++num_distinct_colors;
                        if (num_unitigs > 0) bvb.set(num_unitigs - 1, 1);
                        colors_builder.add_color_set(colors.data, colors.size); // compress colors
                    }
                    bvb.push_back(0);
                    out << ">\n";
                    out.write(unitig.data, unitig.size);
                    out << '\n';
                    num_unitigs += 1;
                } catch (std::exception const& e) {
                    std::cerr << e.what() << std::endl;
                    exit(1);
                }
            }
        );
        out.close();

        if (num_unitigs < std::numeric_limits<uint32_t>::max()) 
            throw std::runtime_error("Number of unitigs > than allowed maximum of " + std::to_string(std::numeric_limits<uint32_t>::max()));
        m_u2c = ranked_bit_vector(std::move(bvb));
        colors_builder.build(m_ccs);

        if (build_parameters.verbose) {
            std::cerr << "num_unitigs " << num_unitigs << "\n";
            std::cerr << "num_distinct_colors " << num_distinct_colors << "\n";
            std::cerr << "m_u2c.size() " << m_u2c.size() << "\n";
            std::cerr << "m_u2c.size1() " << m_u2c.size1() << "\n";
            std::cerr << "m_u2c.size0() " << m_u2c.size0() << "\n";
        }
    }

    { // FIXME
        if (build_parameters.verbose) std::cerr << "step 2. build m_k2u\n";
        sshash::build_configuration sshash_config;
        sshash_config.k = builder.config.k;
        sshash_config.m = builder.config.m;
        sshash_config.canonical_parsing = builder.config.canonical_parsing;
        sshash_config.verbose = builder.config.verbose;
        sshash_config.tmp_dirname = builder.config.tmp_dirname;
        sshash_config.print();
        m_k2u.build(builder.config.file_base_name + ".fa", sshash_config);
        try {  // remove unitig file
            std::remove((builder.config.file_base_name + ".fa").c_str());
        } catch (std::exception const& e) { std::cerr << e.what() << std::endl; }
    }

    { // FIXME
        essentials::logger("step 3. write filenames");
        m_filenames.build(builder.ggcat->filenames());
    }

    if (build_parameters.check) {
        std::cerr << "step 4. check correctness...\n";
        builder.ggcat->loop_through_unitigs(
            [&](ggcat::Slice<char> const unitig, 
                ggcat::Slice<uint32_t> const colors,
                bool) // everything has the same color
            {
                auto lookup_result = m_k2u.lookup_advanced(unitig.data);
                uint32_t unitig_id = lookup_result.contig_id;
                uint32_t color_id = u2c(unitig_id);
                for (std::size_t i = 1; i != unitig.size - m_k2u.k() + 1; ++i) {
                    uint32_t got = m_k2u.lookup_advanced(unitig.data + i).contig_id;
                    if (got != unitig_id) throw std::runtime_error("got unitig_id " + std::to_string(got) + " but expected " + std::to_string(unitig_id));
                }
                auto fwd_it = m_ccs.colors(color_id);
                uint64_t size = fwd_it.size();
                if (size != colors.size) throw std::runtime_error("got colors list of size " + std::to_string(size) + " but expected " + std::to_string(colors.size));
                for (uint64_t i = 0; i != size; ++i, ++fwd_it) {
                    uint32_t ref = *fwd_it;
                    if (ref != colors.data[i]) throw std::runtime_error("got ref " + std::to_string(ref) + " but expected " + std::to_string(colors.data[i]));
                }
            },
            build_parameters.nthreads);
        std::cerr << "DONE!\n";
    }
}

} // namespace kaminari

#endif // KAMINARI_INDEX_HPP