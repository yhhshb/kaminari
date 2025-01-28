#ifndef MINIMIZER_HPP
#define MINIMIZER_HPP

#include <vector>
#include <string>
#include "../bundled/biolib/include/constants.hpp"

namespace minimizer {

#pragma pack(push, 2)
struct record_t {
    uint64_t itself;
    uint64_t id;
    uint32_t p1;
    uint32_t size;
};
#pragma pack(pop)

struct mm_quartet_t {
    mm_quartet_t() : hash(0), itself(0), p1(0), size(0) {}
    void clear() noexcept {
        itself = 0;  // AAA...A
        hash = std::numeric_limits<decltype(hash)>::max();
        p1 = std::numeric_limits<decltype(p1)>::max();
        size = std::numeric_limits<decltype(size)>::max();
    }

    friend std::ostream& operator<<(std::ostream& os, mm_quartet_t const& other) {
        return os << other.itself << " " << other.id << " " << other.p1 << " "
                  << other.size;
    } //removed hash from output because __uint128_t not supported

    __uint128_t hash;  // minimizer hash
    uint64_t id;
    uint64_t itself;  // 2-bit minimizer itself
    uint32_t p1;       // position inside first k-mer of the super-k-mer
    uint32_t size;     // size (number of k-mers) in the super-k-mer
};

template <typename MinimizerHasher, class VectorType>
uint64_t from_string(
    char const* contig, 
    std::size_t contig_size, 
    uint32_t k,
    uint32_t m, 
    uint64_t seed, 
    bool canonical_m_mers,
    std::size_t& mm_count,
    VectorType& accumulator) 
{
    static_assert(std::is_same<typename VectorType::value_type, record_t>::value);
    std::size_t buf_pos, min_pos;
    mm_quartet_t current;
    uint64_t shift = 2 * (m - 1);
    uint64_t mask = (1ULL << (2 * m)) - 1;
    uint64_t mm[2] = {0, 0};
    uint64_t nbases_since_last_break = 0;
    uint32_t sks = 0, p1 = 0;
    uint64_t kmer_count;
    std::vector<mm_quartet_t> buffer(k - m + 1); //minmers / 1kmer
    int c;
    uint8_t z;
    bool find_brand_new_min = false;

    auto update_output = [](decltype(accumulator)& accumulator, mm_quartet_t const& added) {
        accumulator.push_back({added.itself, added.id, added.p1, added.size});
    };

    [[maybe_unused]] auto bit_to_nuc = [](uint64_t bits, uint32_t len) {
        //used for debugging
        std::string nuc;
        for (uint32_t i = 0; i < len; ++i) {
            uint8_t b = (bits >> (2 * (len - i - 1))) & 3;
            switch (b) {
                case 0: nuc += 'A'; break;
                case 1: nuc += 'C'; break;
                case 2: nuc += 'G'; break;
                case 3: nuc += 'T'; break;
            }
        }
        return nuc;
    };

    assert(k >= m);

    buf_pos = 0;
    min_pos = buffer.size();
    kmer_count = 0;
    mm_count = 0;
    z = 0;
    std::array<uint64_t, 2> hash_output = {0, 0};


    for (size_t i = 0; i < contig_size; ++i) {
        //std::cerr << "i " << i << " contig[i] " << contig[i] << "\n";
        c = constants::seq_nt4_table[static_cast<uint8_t>(contig[i])];
        current.clear();
        if (c < 4) [[likely]] { //a t c or g (caps or not)
                mm[0] = (mm[0] << 2 | c) & mask;            // forward m-mer
                mm[1] = (mm[1] >> 2) | (3ULL ^ c) << shift; // reverse m-mer
                if (canonical_m_mers && mm[0] != mm[1]) z = mm[0] < mm[1] ? 0 : 1;  // strand, if symmetric k-mer then use previous strand
                ++nbases_since_last_break;
                if (nbases_since_last_break >= m) { //at least size of minmer (m)
                    /*std::cerr << "\tnb_bases >= m, \n" <<
                                    "\t\t forward : " << mm[0] << " = " << bit_to_nuc(mm[0], m) << "\n" <<
                                    "\t\t reverse : " << mm[1] << " = " << bit_to_nuc(mm[1], m) << "\n" <<
                                    "\t\t chosen : " << mm[z] << " = " << bit_to_nuc(mm[z], m) << "\n";*/
                    current.itself = mm[z]; // itself = canonical encoded
                    
                    //hash_output = MinimizerHasher::hash(mm[z], seed);
                    //current.hash = (static_cast<__uint128_t>(hash_output[1]) << 64) | hash_output[0];
                    current.hash = MinimizerHasher::hash(mm[z], seed);  // insert new hash inside buffer (murmurhash)
                    //std::cerr << "\t\t hash : " << current.hash << "\n";
                    current.p1 = i - m + 1;  // FIXME this is NOT the position inside the super-k-mer!
                    current.id = mm_count++;
                    if (nbases_since_last_break == k) ++kmer_count;
                    if (nbases_since_last_break == k + 1) [[unlikely]] {  
                        // happens once, after seeing all the m-mers of the first kmer, time to check which one is the smallest (the minmer)
                        //std::cerr << "unlikely bases_since_last_break == k + 1\n";
                        min_pos = p1 = 0;
                        //smallest minmer according to murmurhash
                        for (std::size_t j = 0; j < buffer.size(); ++j) {
                            if (buffer[j].hash < buffer[min_pos].hash) {
                                min_pos = j;
                            }
                        }
                        p1 = min_pos;
                        //std::string tmp(contig + p1, m);
                        //std::cerr << "=== Smallest minmer hash = " <<  buffer[p1].hash << " itself : " << buffer[p1].itself << " which corresponds to " << tmp << "\n";
                        sks = 1; // number of k-mers after a break is 1
                    }
                    if (nbases_since_last_break >= k + 1) [[likely]] {  // time to update the minimum, if necessary
                        //std::cerr << "likely bases_since_last_break >= k + 1\n";
                        assert(sks != 0);
                        assert(sks <= k - m + 1); //sks is number of kmer in the superkmer
                        if (((buf_pos) % buffer.size()) == min_pos) {  // old minimum outside window, time to save it and go next
                            buffer[min_pos].p1 = p1;
                            buffer[min_pos].size = sks;
                            //std::cerr << "update output with : " << buffer[min_pos].itself << " = " << bit_to_nuc(buffer[min_pos].itself, m) << "\n";
                            update_output(accumulator, buffer[min_pos]);  // we save the old minimum, length on the right is k by definition
                            sks = 0;
                            find_brand_new_min = true;  // also update p1
                        } else if (current.hash < buffer[min_pos].hash) {
                            buffer[min_pos].p1 = p1;
                            buffer[min_pos].size = sks;
                            //std::cerr << "update output(2) with : " << buffer[min_pos].itself << " = " << bit_to_nuc(buffer[min_pos].itself, m) << "\n";
                            update_output(accumulator, buffer[min_pos]);  // new minimum
                            sks = 0;
                            p1 = k - m;
                            min_pos = buf_pos;  // actual update is outside if
                        }
                        ++sks;
                        ++kmer_count;

                    }
                    buffer[buf_pos++] = current;
                    buf_pos %= buffer.size();  // circular buffer
                    if (find_brand_new_min) {  // find new minimum if the old one dropped out the window
                        find_brand_new_min = false;
                        min_pos = buf_pos;
                        p1 = 0;
                        uint32_t tmp = 1;
                        for (std::size_t j = (buf_pos + 1) % buffer.size(); j < buffer.size(); ++j) {
                            if (buffer[min_pos].hash > buffer[j].hash) {
                                min_pos = j;
                                p1 = tmp;
                            }
                            ++tmp;
                        }
                        for (std::size_t j = 0; j <= buf_pos; ++j) {
                            if (buffer[min_pos].hash > buffer[j].hash) {
                                min_pos = j;
                                p1 = tmp;
                            }
                            ++tmp;
                        }
                    }
                }
            }
        else [[unlikely]] { // N or other character
            nbases_since_last_break = 0;
            if (min_pos < buffer.size()) {
                buffer[min_pos].p1 = p1;
                buffer[min_pos].size = sks;
                //std::cerr << "update output(3) with : " << buffer[min_pos].itself << " = " << bit_to_nuc(buffer[min_pos].itself, m) << "\n";
                update_output(accumulator, buffer[min_pos]);  // push current minimum if available
            }
            sks = 0;  // impossible value, wait for reinitialization of the first window
            min_pos = buffer.size();
            buf_pos = 0; // we always restart at the beginning of the buffer -> this allows to
                         // use min_pos as the position of the minimizer inside the first k-mer
        }
    }
    if (nbases_since_last_break == k) {  // contig.length == 1
        min_pos = p1 = 0;
        sks = 1;
        for (std::size_t j = 0; j < buffer.size(); ++j) {
            if (buffer[j].hash < buffer[min_pos].hash) {
                min_pos = j;
                p1 = min_pos;
            }
        }
    }
    if (min_pos < buffer.size()) {
        buffer[min_pos].p1 = p1;
        buffer[min_pos].size = sks;
        //std::cerr << "update output(4) with : " << buffer[min_pos].itself << " = " << bit_to_nuc(buffer[min_pos].itself, m) << "\n";
        update_output(accumulator, buffer[min_pos]);  // push last minimum if available
        sks = 1;
    }
    return kmer_count;
}

}  // namespace minimizer

#endif // MINIMIZER_HPP