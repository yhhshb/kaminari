#ifndef IO_HELPER_HPP
#define IO_HELPER_HPP

#include "../../bundled/Minimizers/src/files/stream_reader_library.hpp"
#include "../../bundled/Minimizers/src/files/stream_writer_library.hpp"

#include <vector>
#include <string>
#include <memory>
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <filesystem> 

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
        // FIX: Check existence BEFORE calling the strict library allocator
        // The library calls exit(1) if the file is missing, so we must protect it.
        if (!std::filesystem::exists(filename)) {
            // File does not exist -> Leave _reader as null, set EOF true.
            // is_open() will return false, handling this gracefully.
            _reader.reset(nullptr);
            _eof_reached = true;
            return;
        }

        // Allocate reader using the factory (handles .lz4 extension automatically)
        _reader.reset(stream_reader_library::allocate(filename));
        
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
        if (!ensure_available(2)) return false;

        // 1. Read Minimizer
        out_minimizer = _buffer[_pos];

        // 2. Peek Header to calculate size
        uint64_t header = _buffer[_pos + 1];
        
        // Calculate total payload words (including the header word itself)
        uint64_t colors_per_word = 64 / bits_per_color;
        uint64_t list_size = header >> (64 - bits_per_color);
        size_t payload_words = (list_size + colors_per_word) / colors_per_word;
        
        // Total words to advance = 1 (Minimizer) + payload_words
        size_t total_len = 1 + payload_words;

        // 3. Advance logic
        if (_pos + total_len <= _limit) {
            _pos += total_len;
        } else {
            size_t available = _limit - _pos;
            size_t remaining_to_skip = total_len - available;
            
            _pos = _limit; // Buffer fully consumed

            while (remaining_to_skip > 0) {
                if (_pos >= _limit) {
                    if (_eof_reached) return false; 
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
     */
    bool next_element(uint64_t& out_minimizer, std::vector<uint64_t>& out_payload, uint64_t bits_per_color) {
        if (!ensure_available(2)) return false;

        out_minimizer = _buffer[_pos++]; 
        uint64_t header = _buffer[_pos];
        uint64_t colors_per_word = 64 / bits_per_color;
        uint64_t list_size = header >> (64 - bits_per_color);
        size_t payload_words = (list_size + colors_per_word) / colors_per_word;

        if (!ensure_available(payload_words)) return false; 

        out_payload.resize(payload_words);
        std::memcpy(out_payload.data(), &_buffer[_pos], payload_words * sizeof(uint64_t));
        _pos += payload_words;

        return true;
    }

private:
    std::unique_ptr<stream_reader> _reader;
    std::vector<uint64_t> _buffer;
    size_t _pos;    
    size_t _limit;  
    bool _eof_reached;

    bool ensure_available(size_t n) {
        if (_pos + n <= _limit) return true;
        if (_eof_reached) return false;

        size_t remaining = _limit - _pos;
        if (remaining > 0) {
            std::memmove(_buffer.data(), &_buffer[_pos], remaining * sizeof(uint64_t));
        }
        _pos = 0;
        _limit = remaining;

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
        ensure_available(1); 
    }
};

class BinaryFileIterator {
public:
    using iterator_category = std::input_iterator_tag;
    using value_type        = uint64_t;
    using difference_type   = std::ptrdiff_t;
    using pointer           = const uint64_t*;
    using reference         = const uint64_t&;

    // Constructor
    BinaryFileIterator(const std::string& filename, size_t buffer_size_bytes = 4 * 1024 * 1024)
        : _buffer_pos(0), _is_end(false) 
    {
        // Allocate the stream on the heap and manage via shared_ptr
        _stream = std::make_shared<std::ifstream>(filename, std::ios::binary);
        
        if (!_stream->is_open()) {
            _is_end = true;
        } else {
            _buffer.resize(buffer_size_bytes / sizeof(uint64_t));
            fill_buffer();
        }
    }

    // Default constructor (End Iterator)
    BinaryFileIterator() 
        : _buffer_pos(0), _is_end(true) 
    {}

    uint64_t operator*() const {
        return _buffer[_buffer_pos];
    }

    BinaryFileIterator& operator++() {
        _buffer_pos++;
        if (_buffer_pos >= _items_in_buffer) {
            fill_buffer();
        }
        return *this;
    }

    bool operator!=(const BinaryFileIterator& other) const {
        return _is_end != other._is_end;
    }
    
    BinaryFileIterator(const BinaryFileIterator&) = default;
    BinaryFileIterator& operator=(const BinaryFileIterator&) = default;

private:
    std::shared_ptr<std::ifstream> _stream; 
    std::vector<uint64_t> _buffer;
    size_t _buffer_pos;
    size_t _items_in_buffer = 0;
    bool _is_end;

    void fill_buffer() {
        if (!_stream || !_stream->is_open()) { 
            _is_end = true;
            return;
        }
        
        _stream->read(reinterpret_cast<char*>(_buffer.data()), _buffer.size() * sizeof(uint64_t));
        std::streamsize bytes_read = _stream->gcount(); 

        if (bytes_read == 0) {
            _is_end = true;
            return;
        }

        _items_in_buffer = bytes_read / sizeof(uint64_t);
        _buffer_pos = 0;
    }
};

} // namespace helpers
} // namespace kaminari

#endif // IO_HELPER_HPP