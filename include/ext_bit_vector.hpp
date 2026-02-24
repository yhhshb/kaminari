#ifndef EXT_BIT_VECTOR_HPP
#define EXT_BIT_VECTOR_HPP

#include <vector>
#include <string>
#include <cstdint>
#include <fcntl.h>
#include <unistd.h>
#include <system_error>
#include <cerrno>

#include <iostream>


//LSB-first packing inside each uint64_t
class ext_bit_vector {
public:
    static constexpr size_t BITS = 64;

    ext_bit_vector(){};

    ext_bit_vector(std::string basename, size_t ram_limit_bytes = 256ull<<20)
        : m_ram_limit_bytes(ram_limit_bytes),
          m_size(0),
          m_fd(nullptr),
          m_flushed_words(0)
    {
        m_filename = basename + ".data";
        open_file();
    }

    // number of bits
    size_t size() const { return m_size; }

    // resize (new bits are zero)
    //considering m_ram_limit_bytes is enough to handle 1 colorset
    void expand(size_t n_bits) {
        uint64_t last = m_data.back();
        if (m_data.size() * sizeof(uint64_t) + (n_bits+7)/8 >= m_ram_limit_bytes) {
            flush_available((m_size % BITS != 0));
            if (m_size % BITS != 0) {
                m_data.push_back(last); // keep last partial word
            }
        }
        m_size += n_bits;
        size_t n_RAM_words = (m_size + BITS - 1) / BITS - m_flushed_words;
        if (n_RAM_words > m_data.size()) {
            m_data.resize(n_RAM_words, 0);
        }
    }

    // append one bit
    void push_back(bool bit) {
        size_t idx = (m_size / BITS) - m_flushed_words;
        size_t pos = m_size % BITS;
        if (pos == 0){
            if (m_data.size() * sizeof(uint64_t) + sizeof(uint64_t) >= m_ram_limit_bytes) {
                flush_available();
            }
            m_data.push_back(0);
        } 
        if (bit) m_data[idx] |= (uint64_t(1) << pos);
        ++m_size;
    }


    // append multiple bits contained inside data (lower n_bits are used)
    void push_back(uint64_t data, uint64_t n_bits) {
        if (n_bits == 0) return; 

        size_t pos = m_size % BITS;
        
        // Handle Auto-Growth for new words
        if (pos == 0) {
            if (m_data.size() * sizeof(uint64_t) + sizeof(uint64_t) >= m_ram_limit_bytes) {
                flush_available();
            }
            m_data.push_back(0);
        } 

        size_t idx = (m_size / BITS) - m_flushed_words;

        // Write the lower part of data
        m_data[idx] |= data << pos; 
        m_size += n_bits;

        if (pos + n_bits > BITS) {
            uint64_t remaining_bits = data >> (BITS - pos);

            // If the current word (idx) was the last one, and we need to add a new one,
            // check RAM limit first.
            if (m_data.size() * sizeof(uint64_t) + sizeof(uint64_t) >= m_ram_limit_bytes) {
                // Flush but KEEP the word we just wrote the first part into
                flush_available(true); 
                // After flush(true), m_data has 1 element (our partial word).
                // The next part goes into the *new* word at index 1.
                idx = 0; 
            }

            if (m_data.size() <= idx + 1) {
                m_data.push_back(0);
            }

            m_data[idx + 1] |= remaining_bits;
        }
    }

    // set bit
    void set(size_t i, bool bit) {
        if (i >= m_size) throw std::out_of_range("ext_bit_vector::set");
        size_t idx = (i / BITS) - m_flushed_words; //here danger, idx has to be >= 0 but should be the case as set() is only used after resize for semi-dense colorset
        size_t pos = i % BITS;

        if (bit) m_data[idx] |= (uint64_t(1) << pos);
        else     m_data[idx] &= ~(uint64_t(1) << pos);
    }

    // get bit, obsolete
    bool get(size_t i) const {
        if (i >= m_size) throw std::out_of_range("ext_bit_vector::get");
        size_t idx = (i / BITS) - m_flushed_words;
        size_t pos = i % BITS;
        return (m_data[idx] >> pos) & 1;
    }

    // flush in-memory data to disk
    //TODO : add to build()
    void flush_all() {
        if (m_fd == nullptr) open_file();
        if (m_data.empty()) return;
        size_t written = fwrite(m_data.data(), sizeof(uint64_t), m_data.size(), m_fd);
        if (written != m_data.size()) {
            fclose(m_fd);
            m_fd = nullptr;
            throw std::runtime_error("Error writing to file " + m_filename);
        }
        m_flushed_words += written;
        m_data.clear();
    }


private:
    size_t m_ram_limit_bytes;
    size_t m_size; // number of valid bits
    FILE* m_fd; //file descriptor
    size_t m_flushed_words; // number of words already flushed to disk
    std::string m_filename;
    std::vector<uint64_t> m_data;  // packed storage


    void open_file() {
        if (m_fd != nullptr) return;
        m_fd = fopen(m_filename.c_str(), "a+b");
        if (m_fd == nullptr) {
            throw std::runtime_error("Cannot open file " + m_filename);
        }
    }

    void flush_available(bool keep_last = false) {
        if (m_fd == nullptr) open_file();
        
        if (m_data.empty()) return;

        size_t to_flush = m_data.size();
        if (keep_last) {
            if (to_flush == 0) return; // Should not happen but safe
            to_flush--;
        }

        if (to_flush > 0) {
            size_t written = fwrite(m_data.data(), sizeof(uint64_t), to_flush, m_fd);
            if (written != to_flush) {
                fclose(m_fd);
                m_fd = nullptr;
                throw std::runtime_error("Error writing to file " + m_filename);
            }
            m_flushed_words += written;
        }

        // Critical: Handle the kept word
        if (keep_last) {
            uint64_t last = m_data.back();
            m_data.clear();
            m_data.push_back(last);
        } else {
            m_data.clear();
        }
    }
    
};

#endif // EXT_BIT_VECTOR_HPP