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
            uint64_t ram_limit_bytes = 256ull<<20)
        : m_basename(basename),
          m_total_elements(total_elements),
          m_width(width),
          m_mask((width == 64) ? uint64_t(-1) : ((uint64_t(1) << width) - 1)),
          m_first_chunk_index(0)
        {
            if (m_width == 0 || m_width > 64)
                throw std::runtime_error("Invalid width");

            // how many words in buffer
            uint64_t words_needed = words_for(m_total_elements * m_width);
            uint64_t words_fit = ram_limit_bytes / sizeof(uint64_t);
            uint64_t words_init = std::min(words_needed, words_fit);
            m_data.resize(words_init, 0);
            m_buffer_bits    = words_init * 64;
            m_chunk_capacity = m_buffer_bits / m_width;

            // open data file (C IO)
            m_out = std::fopen((m_basename + ".data").c_str(), "wb");
            if (!m_out) throw std::runtime_error("Failed to open .data file");
        }

        // initialize a new chunk, returns how many elements to insert
        uint64_t init_chunk() {
            m_filled = 0;
            // how many elements needed to write before full chunk 
            return std::min(m_chunk_capacity, m_total_elements);
        }

        void set(uint64_t idx, uint64_t v) {
            assert(idx >= m_first_chunk_index &&
                idx < m_first_chunk_index + m_chunk_capacity);

            uint64_t rel_idx = idx - m_first_chunk_index;
            uint64_t pos     = rel_idx * m_width;
            uint64_t block   = pos >> 6;
            uint64_t shift   = pos & 63;

            m_data[block] &= ~(m_mask << shift);
            m_data[block] |= (v & m_mask) << shift;

            uint64_t res_shift = 64 - shift;
            if (res_shift < m_width) {
                m_data[block + 1] &= ~(m_mask >> res_shift);
                m_data[block + 1] |= (v & m_mask) >> res_shift;
            }
            m_filled++;
        }


        uint64_t flush() {
            uint64_t last_word_index = (m_filled * m_width) % 64;
            uint64_t used_words = ((m_filled - m_first_chunk_index) * m_width + 63) / 64; 

            if (last_word_index != 0) { //overlap in chunk last word
                if (used_words == 0){ //vector final word, not fully used
                    m_written += std::fwrite(m_data.data(), sizeof(uint64_t), 1, m_out);
                    
                } else {
                    // Write all but last word
                    m_written += std::fwrite(m_data.data(), sizeof(uint64_t), used_words - 1, m_out);
                    // Carry last word
                    uint64_t last_word = m_data[used_words - 1];
                    m_data[0] = last_word;
                }
                m_chunk_capacity = (m_buffer_bits - (64-last_word_index)) / m_width;
                
            } else { //no overlap
                // Write everything
                m_written += std::fwrite(m_data.data(), sizeof(uint64_t), used_words, m_out);
                m_chunk_capacity = m_buffer_bits / m_width;
            }

            m_first_chunk_index = m_filled;
            return std::min(m_chunk_capacity, m_total_elements - m_first_chunk_index);
        }

        void build() {
            flush();
            if (m_out) std::fclose(m_out);

            // write metadata (C IO)
            header_t header{m_total_elements, m_width, m_mask};
            FILE* meta = std::fopen((m_basename + ".metadata").c_str(), "wb");
            if (!meta) throw std::runtime_error("Failed to open metadata file");
            std::fwrite(&header, sizeof(header), 1, meta);
            std::fclose(meta);

            std::cerr << 
                "ColorMapper stats:\n" <<
                "\tcontains " << m_total_elements << " slots (1/minimizer)\n" <<
                "\tbits/slot = " << m_width << "\n" <<
                "\ttotal bits = " << words_for(m_total_elements * m_width)*64 << " (" << words_for(m_total_elements * m_width)*64/8.0/1048576.0 << " MB)\n";
        }

        uint64_t chunk_capacity() const { return m_chunk_capacity; }
        uint64_t total_size() const { return m_total_elements; }

    private:
        std::string m_basename;
        uint64_t m_total_elements;
        uint64_t m_width;
        uint64_t m_mask;

        uint64_t m_buffer_bits;
        uint64_t m_chunk_capacity;
        uint64_t m_first_chunk_index;
        uint64_t m_filled = 0;
        uint64_t m_written = 0;

        std::vector<uint64_t> m_data;
        FILE* m_out = nullptr;
    };


    // -------- COMPACT VECTOR ----------
    colormapper() : m_total_elements(0), m_width(0), m_mask(0) {}

    static colormapper load(const std::string& basename) {
        colormapper cv;

        // read metadata (C IO)
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