#ifndef COLORSETS
#define COLORSETS

#include "ext_bit_vector.hpp"
#include "codes.hpp"

#include "../bundled/biolib/include/elias_fano.hpp"
#include "../bundled/biolib/include/io.hpp"

struct colorsets {
    struct builder {
        builder(std::string basename, size_t num_docs, size_t ram_limit_bytes = 256ull<<20)
            : m_num_docs(num_docs),
            m_sparse_set_threshold_size(0.25 * num_docs),
            m_dense_set_threshold_size(0.75 * num_docs),
            m_offsets(),
            m_basename(std::move(basename)),
            m_colors(m_basename, ram_limit_bytes)
        {
            m_building_offsets.push_back(0);

            std::cerr << "Creating colorsets_builder with " << num_docs << " documents\n" << 
                            "sparse set threshold size: " << m_sparse_set_threshold_size << "\n" <<
                            "dense set threshold size: " << m_dense_set_threshold_size << "\n";
        }

        void add_color_set(uint32_t const * const colors, size_t list_size) {
            // encode list_size
            delta_encoder(m_colors, list_size);

            if (list_size < m_sparse_set_threshold_size) {
                // SPARSE SET: gaps + delta
                uint32_t prev_val = colors[0];
                delta_encoder(m_colors, prev_val);
                for (size_t i = 1; i != list_size; ++i) {
                    uint32_t val = colors[i];
                    if (val < prev_val + 1) throw std::runtime_error("colorsets_builder::add_color_set: colors not sorted or repeated");
                    delta_encoder(m_colors, val - (prev_val + 1));
                    prev_val = val;
                }
            } else if (list_size < m_dense_set_threshold_size) {
                // MIXED SET: bitmap
                size_t old_size = m_colors.size();
                m_colors.expand(m_num_docs);
                for (size_t i = 0; i != list_size; ++i) {
                    m_colors.set(old_size + colors[i], true);
                }
            } else {
                // DENSE SET: complementary gaps + delta
                bool first = true;
                uint32_t val = 0;
                uint32_t prev_val = -1;
                uint32_t written = 0;
                for (size_t i = 0; i != list_size; ++i) {
                    uint32_t x = colors[i];
                    while (val < x) {
                        if (first) {
                            delta_encoder(m_colors, val);
                            first = false;
                            ++written;
                        } else {
                            if (val < prev_val + 1) throw std::runtime_error("colorsets_builder::add_color_set: colors not sorted or repeated");
                            delta_encoder(m_colors, val - (prev_val + 1));
                            ++written;
                        }
                        prev_val = val;
                        ++val;
                    }
                    val++;  // skip x
                }
                while (val < m_num_docs) {
                    if (val < prev_val + 1) throw std::runtime_error("colorsets_builder::add_color_set: colors not sorted or repeated");
                    delta_encoder(m_colors, val - (prev_val + 1));
                    prev_val = val;
                    ++val;
                    ++written;
                }
                if (written != m_num_docs - list_size) 
                    throw std::runtime_error("colorsets_builder::add_color_set: unexpected number of written values");
            }
            m_building_offsets.push_back(m_colors.size());  
        }

        void build() {
            m_colors.flush_all();
                        
            m_offsets = bit::ef::array(
                m_building_offsets.begin(), 
                m_building_offsets.size(),
                m_building_offsets.back()); //max possible value
            write_metadata();

            //free memory
            m_building_offsets.clear();
            m_building_offsets.shrink_to_fit();

            std::cerr << 
                "Colorsets_builder stats:\n" <<
                "\tprocessed " << m_offsets.size() - 1 << " unique color sets\n" <<
                "\tstored " << m_colors.size() << " bits (" << (m_colors.size()/8.0)/1048576.0 << " MB)\n" <<
                "\tstored " << m_offsets.bit_size() << " bits for offsets (" << (m_offsets.bit_size()/8.0)/1048576.0 << " MB)\n" <<
                "\ttotal bits = " << (m_offsets.bit_size() + m_colors.size()) << " (" << (m_offsets.bit_size() + m_colors.size())/8.0/1048576.0 << " MB)\n";
        }

        private:
            void write_metadata() const {
                std::ofstream metadata_out(m_basename + ".metadata", std::ios::binary);
                if (!metadata_out) throw std::runtime_error("Failed to open metadata file for writing");

                io::saver visitor(metadata_out);
                visitor.visit(m_num_docs);
                visitor.visit(m_sparse_set_threshold_size);
                visitor.visit(m_dense_set_threshold_size);
                visitor.visit(m_offsets);
            }

            size_t m_num_docs;
            size_t m_sparse_set_threshold_size;
            size_t m_dense_set_threshold_size;
            bit::ef::array m_offsets;
            std::string m_basename;
            ext_bit_vector m_colors;
            std::vector<size_t> m_building_offsets;
            ext_bit_parser m_parser;
    }; //end builder

    // ---------- COLORSETS ----------------------
    colorsets() : 
        m_num_docs(0), 
        m_sparse_set_threshold_size(0), 
        m_dense_set_threshold_size(0) 
    {}

    static void load(colorsets& cs, const std::string& basename) {
        std::ifstream metadata_in(basename + ".metadata", std::ios::binary);
        if (!metadata_in) throw std::runtime_error("Failed to open metadata file for reading");

        io::loader visitor(metadata_in);
        visitor.visit(cs.m_num_docs);
        visitor.visit(cs.m_sparse_set_threshold_size);
        visitor.visit(cs.m_dense_set_threshold_size);
        visitor.visit(cs.m_offsets);

        cs.m_mapped_file = mymm::immap<uint64_t>(basename + ".data", 0, 0);
    }


    uint64_t m_num_docs;
    uint64_t m_sparse_set_threshold_size;
    uint64_t m_dense_set_threshold_size;
    bit::ef::array m_offsets;
    mymm::immap<uint64_t> m_mapped_file;
};


#endif // COLORSETS