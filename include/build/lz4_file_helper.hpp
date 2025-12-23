#ifndef LZ4_FILE_HELPER_HPP
#define LZ4_FILE_HELPER_HPP

#include "../../bundled/Minimizers/src/files/stream_reader_library.hpp"
#include "../../bundled/Minimizers/src/files/stream_writer_library.hpp"

#include <vector>
#include <string>
#include <memory>
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <iostream>

namespace kaminari {
namespace helpers {

class SparseFileParser {
public:
    /**
     * @brief Construct a new Sparse File Parser object
     * * @param filename Path to the .lz4 sparse file
     * @param buffer_size_mb Internal buffer size in MB (default 4MB)
     */
    SparseFileParser(const std::string& filename, size_t buffer_size_mb = 4) 
        : _pos(0), _limit(0), _eof_reached(false) 
    {
        // Allocate reader using the factory (handles .lz4 extension automatically)
        _reader.reset(stream_reader_library::allocate(filename));
        
        // Handle case where file cannot be opened (e.g. sparse file doesn't exist)
        if (!_reader || !_reader->is_open()) {
            _eof_reached = true; 
        } else {
            // Allocate buffer: ensure it's a multiple of uint64_t and at least 1KB
            size_t num_words = (buffer_size_mb * 1024 * 1024) / sizeof(uint64_t);
            _buffer.resize(std::max<size_t>(1024, num_words));
            refill();
        }
    }

    // Delete copy constructors to prevent stream state duplication issues
    SparseFileParser(const SparseFileParser&) = delete;
    SparseFileParser& operator=(const SparseFileParser&) = delete;

    /**
     * @brief Checks if the underlying file was successfully opened.
     */
    bool is_open() const { 
        return _reader && _reader->is_open(); 
    }

    /**
     * @brief Step 1 Helper: Reads the next minimizer and efficiently skips its payload.
     * * @param out_minimizer Reference to store the read minimizer.
     * @param bits_per_color Used to calculate payload size from the header.
     * @return true if a minimizer was read, false if EOF.
     */
    bool next_minimizer_only(uint64_t& out_minimizer, uint64_t bits_per_color) {
        // We need at least 2 words to make a decision: Minimizer + Header
        if (!ensure_available(2)) return false;

        // 1. Read Minimizer
        out_minimizer = _buffer[_pos];

        // 2. Peek Header to calculate size
        uint64_t header = _buffer[_pos + 1];
        
        // Calculate total payload words (including the header word itself)
        // Formula: (list_size + capacity_per_word) / capacity_per_word
        uint64_t colors_per_word = 64 / bits_per_color;
        uint64_t list_size = header >> (64 - bits_per_color);
        size_t payload_words = (list_size + colors_per_word) / colors_per_word;
        
        // Total words to advance = 1 (Minimizer) + payload_words
        size_t total_len = 1 + payload_words;

        // 3. Advance logic
        // If the entire element fits in the current buffer, just jump.
        if (_pos + total_len <= _limit) {
            _pos += total_len;
        } else {
            // The element is split across the buffer boundary.
            // Since we only want to SKIP, we don't need to contiguous-load the payload.
            // Just consume what is here and skip the rest from stream/refill.
            
            // Consume whatever is currently in buffer
            size_t available = _limit - _pos;
            size_t remaining_to_skip = total_len - available;
            
            _pos = _limit; // Buffer fully consumed

            // Skip remaining words (handling refills if payload > buffer, though unlikely)
            while (remaining_to_skip > 0) {
                if (_pos >= _limit) {
                    if (_eof_reached) return false; // Unexpected EOF inside element
                    refill();
                    if (_limit == 0) return false;
                }
                size_t skip_now = std::min(remaining_to_skip, _limit - _pos);
                _pos += skip_now;
                remaining_to_skip -= skip_now;
            }
        }

        return true;
    }

    /**
     * @brief Step 4 Helper: Reads the next minimizer and extracts the full payload.
     * * @param out_minimizer Reference to store the read minimizer.
     * @param out_payload Vector to store the payload (will be resized).
     * @param bits_per_color Used to calculate payload size.
     * @return true if element read, false if EOF.
     */
    bool next_element(uint64_t& out_minimizer, std::vector<uint64_t>& out_payload, uint64_t bits_per_color) {
        // Need Minimizer + Header
        if (!ensure_available(2)) return false;

        // 1. Read Minimizer
        out_minimizer = _buffer[_pos++]; 

        // 2. Read Header
        uint64_t header = _buffer[_pos];
        uint64_t colors_per_word = 64 / bits_per_color;
        uint64_t list_size = header >> (64 - bits_per_color);
        size_t payload_words = (list_size + colors_per_word) / colors_per_word;

        // 3. Ensure full payload is available continuously in memory
        if (!ensure_available(payload_words)) {
            // This happens if the file ends prematurely inside a payload
            return false; 
        }

        // 4. Copy payload
        out_payload.resize(payload_words);
        std::memcpy(out_payload.data(), &_buffer[_pos], payload_words * sizeof(uint64_t));
        _pos += payload_words;

        return true;
    }

private:
    std::unique_ptr<stream_reader> _reader;
    std::vector<uint64_t> _buffer;
    size_t _pos;    // Current cursor
    size_t _limit;  // Valid data end
    bool _eof_reached;

    /**
     * @brief Ensures `n` words are available in the buffer starting at `_pos`.
     * If not, moves remaining data to front and reads more.
     */
    bool ensure_available(size_t n) {
        if (_pos + n <= _limit) return true;
        if (_eof_reached) return false;

        // Move remaining data to start of buffer
        size_t remaining = _limit - _pos;
        if (remaining > 0) {
            std::memmove(_buffer.data(), &_buffer[_pos], remaining * sizeof(uint64_t));
        }
        _pos = 0;
        _limit = remaining;

        // Fill rest of buffer
        size_t cap = _buffer.size() - _limit;
        int read_count = _reader->read((void*)(_buffer.data() + _limit), sizeof(uint64_t), cap);
        
        if (read_count > 0) _limit += read_count;
        if (read_count < (int)cap || _reader->is_eof()) _eof_reached = true;

        return _limit >= n;
    }

    void refill() {
        _pos = 0;
        _limit = 0;
        _eof_reached = false;
        ensure_available(1); // Trigger read
    }
};

} // namespace helpers
} // namespace kaminari

#endif // LZ4_FILE_HELPER_HPP