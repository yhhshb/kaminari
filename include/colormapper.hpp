#ifndef KAMINARI_COLORMAPPER_HPP
#define KAMINARI_COLORMAPPER_HPP

#include <vector>
#include <cmath>
#include <iterator>
#include <cstdint>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <string>
#include <iostream>
#include <fstream>
#include <cstring>
#include "mmap.hpp"
#include <cstdio> // For C IO

template <typename WordType = uint64_t>
static uint64_t words_for(uint64_t bits) {
    uint64_t word_bits = sizeof(WordType) * 8;
    return (bits + word_bits - 1) / word_bits;
}

struct colormapper {
    struct header_t {
        uint64_t m_total_elements;
        uint64_t m_width;
        uint64_t m_mask;
    };

    // -------- STREAMING BUILDER ----------
    struct builder {
        builder(const std::string& basename,
                uint64_t total_elements,
                uint64_t width,
                uint64_t ram_limit_bytes = 256ull<<20
        )
            : m_basename(basename),
              m_total_elements(total_elements),
              m_width(width),
              m_mask((width == 64) ? uint64_t(-1) : ((uint64_t(1) << width) - 1)),
              m_first_chunk_index(0),
              m_written_words(0) // Track words written to disk
        {
            if (m_width == 0 || m_width > 64)
                throw std::runtime_error("Invalid width");

            uint64_t words_fit = ram_limit_bytes / sizeof(uint64_t);
            // Ensure we have at least a few words
            if (words_fit < 2) words_fit = 2; 
            
            m_data.resize(words_fit, 0);
            
            m_out = std::fopen((m_basename + ".data").c_str(), "wb");
            if (!m_out) throw std::runtime_error("Failed to open .data file");
        }

        // Initialize a new chunk
        uint64_t init_chunk() {
            m_filled_in_chunk = 0;
            
            // Calculate how many elements fit in the REMAINING space of the buffer
            // Global bit offset where the next element starts
            uint64_t start_bit = m_first_chunk_index * m_width;
            uint64_t start_word_idx = start_bit / 64;
            uint64_t word_offset = start_word_idx - m_written_words;
            
            // Bits available in the buffer from this offset
            uint64_t bits_available = (m_data.size() - word_offset) * 64 - (start_bit % 64);
            uint64_t capacity = bits_available / m_width;

            return std::min(capacity, m_total_elements - m_first_chunk_index);
        }

        void set(uint64_t idx, uint64_t v) {
            uint64_t global_bit_pos = idx * m_width;
            uint64_t global_word_idx = global_bit_pos / 64;
            
            // Map global word index to local buffer index
            // The buffer starts at global word 'm_written_words'
            uint64_t block = global_word_idx - m_written_words;
            uint64_t shift = global_bit_pos % 64;

            if (block >= m_data.size()) {
                throw std::runtime_error("Colormapper builder buffer overflow");
            }

            m_data[block] &= ~(m_mask << shift);
            m_data[block] |= (v & m_mask) << shift;

            uint64_t res_shift = 64 - shift;
            if (res_shift < m_width) {
                if (block + 1 < m_data.size()) {
                    m_data[block + 1] &= ~(m_mask >> res_shift);
                    m_data[block + 1] |= (v & m_mask) >> res_shift;
                }
            }
            m_filled_in_chunk++;
        }

        uint64_t flush() {
            uint64_t total_processed = m_first_chunk_index + m_filled_in_chunk;
            uint64_t global_bit_end = total_processed * m_width;
            uint64_t global_word_end = (global_bit_end + 63) / 64; // Ceil
            uint64_t words_touched = global_word_end - m_written_words;
            
            bool partial_last = (global_bit_end % 64 != 0);
            uint64_t words_to_write = words_touched;
            
            if (partial_last && words_to_write > 0) {
                words_to_write--; 
            }
            
            if (words_to_write > 0) {
                std::fwrite(m_data.data(), sizeof(uint64_t), words_to_write, m_out);
                
                if (partial_last) {
                    m_data[0] = m_data[words_to_write];
                } else {
                    m_data[0] = 0; // Clean start
                }
                
                size_t fill_start = partial_last ? 1 : 0;
                std::fill(m_data.begin() + fill_start, m_data.end(), 0);
                
                m_written_words += words_to_write;
            }
            
            m_first_chunk_index = total_processed;
            
            return m_filled_in_chunk; 
        }

        void build(size_t verbose) {
            flush(); 
            
            uint64_t global_bit_end = m_first_chunk_index * m_width;
            if (global_bit_end % 64 != 0) {
                // The loop in flush didn't write the last word because it was partial.
                // We must write it now.
                std::fwrite(m_data.data(), sizeof(uint64_t), 1, m_out);
            }

            if (m_out) std::fclose(m_out);

            // write metadata
            header_t header{m_total_elements, m_width, m_mask};
            FILE* meta = std::fopen((m_basename + ".metadata").c_str(), "wb");
            if (!meta) throw std::runtime_error("Failed to open metadata file");
            std::fwrite(&header, sizeof(header), 1, meta);
            std::fclose(meta);

            if (verbose >= 1)
            std::cout << 
                "[II] ColorMapper stats:\n" <<
                "\t\tcontains " << m_total_elements << " slots\n" <<
                "\t\tbits/slot = " << m_width << "\n";
        }

    private:
        std::string m_basename;
        uint64_t m_total_elements;
        uint64_t m_width;
        uint64_t m_mask;

        uint64_t m_first_chunk_index; // Global index of the first element in the current chunk
        uint64_t m_filled_in_chunk = 0; // How many added since init_chunk
        uint64_t m_written_words = 0;   // How many uint64_t words fully written to disk

        std::vector<uint64_t> m_data;
        FILE* m_out = nullptr;
    };


    // -------- COLORMAPPER  ----------
    colormapper() : m_total_elements(0), m_width(0), m_mask(0) {}

    static colormapper load(const std::string& basename) {
        colormapper cv;

        // read metadata 
        header_t header;
        FILE* meta_in = std::fopen((basename + ".metadata").c_str(), "rb");
        if (!meta_in) throw std::runtime_error("Cannot open metadata file");
        size_t read_count = std::fread(&header, sizeof(header), 1, meta_in);
        std::fclose(meta_in);
        if (read_count != 1) throw std::runtime_error("Failed to read metadata");

        cv.m_total_elements = header.m_total_elements;
        cv.m_width   = header.m_width;
        cv.m_mask    = header.m_mask;

        uint64_t words = words_for(cv.m_total_elements * cv.m_width) + 1;
        auto mmap = mymm::immap<uint64_t>(basename + ".data", words, 0);
        cv.m_data = mmap.data();

        return cv;
    }

    uint64_t operator[](uint64_t i) const {
        assert(i < m_total_elements);
        uint64_t pos   = i * m_width;
        uint64_t block = pos / 64;
        uint64_t shift = pos % 64;

        return shift + m_width <= 64
            ? (m_data[block] >> shift) & m_mask
            : ((m_data[block] >> shift) |
               ((m_data[block + 1] << (64 - shift)) & m_mask));
    }

    uint64_t size() const { return m_total_elements; }

    uint64_t num_bytes() const {
        return sizeof(header_t) + m_total_elements * sizeof(uint64_t);
    }

private:
    uint64_t m_total_elements;
    uint64_t m_width;
    uint64_t m_mask;
    uint64_t const* m_data;
};

#endif // KAMINARI_COLORMAPPER_HPP