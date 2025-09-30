#ifndef MINIMIZER_HPP
#define MINIMIZER_HPP

#include <vector>
#include <string>
#include "../../bundled/biolib/include/constants.hpp"

namespace minimizer {

const std::array<uint8_t, 256> seq_nt5_table = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

struct record_t {
    uint64_t hash;
    uint32_t size;
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
    int n_minmer = 0;
    std::vector<uint64_t> liste_mini(contig_size-k+1);

    const int z = k - m;

    uint64_t buffer[z + 1]; // nombre de m-mer/k-mer
    for(int x = 0; x <= z; x += 1)
        buffer[x] = UINT64_MAX;

    uint64_t mask = (1ULL << (2 * m)) - 1; 

    // PREMIER M-MER
    uint64_t current_mmer = 0;
    uint64_t cur_inv_mmer = 0;
    uint64_t canon;
    size_t cnt = 0;

    for(uint32_t x = 0; x < m - 1; x += 1)
    {
        // Encode les nucleotides du premier s-mer
        const uint64_t encoded = ((contig[cnt] >> 1) & 0b11);
        current_mmer <<= 2;
        current_mmer |= encoded; // conversion ASCII => 2bits (Yoann)
        current_mmer &= mask;
        cur_inv_mmer >>= 2;
        cur_inv_mmer |= ( (0x2 ^ encoded) << (2 * (m - 1))); // cf Yoann
        cnt          += 1;
    }


    // PREMIER K-MER
    uint64_t minv = UINT64_MAX;
    for(int m_pos = 0; m_pos <= z; m_pos += 1)
    {
        // On calcule les m-mers
        const uint64_t encoded = ((contig[cnt] >> 1) & 0b11); // conversion ASCII => 2bits (Yoann)
        current_mmer <<= 2;
        current_mmer |= encoded;
        current_mmer &= mask;
        cur_inv_mmer >>= 2;
        cur_inv_mmer |= ( (0x2 ^ encoded) << (2 * (m - 1))); // cf Yoann
        
        if (canonical_m_mers){
            canon  = (current_mmer < cur_inv_mmer) ? current_mmer : cur_inv_mmer;
        } else {
            canon  = current_mmer;
        }

        const uint64_t s_hash = MinimizerHasher::hash(canon, seed);
        buffer[m_pos]         = s_hash; // on memorise le hash du mmer
        minv                  = (s_hash < minv) ? s_hash : minv;
        cnt                  += 1; // on avance dans la sequence de nucleotides
    }

    accumulator.push_back({minv, 1});
    n_minmer += 1;


    //ALL OTHER KMERS
    while (cnt < contig_size) {
        const uint64_t encoded = ((contig[cnt] >> 1) & 0b11);
        current_mmer <<= 2;
        current_mmer |= encoded;
        current_mmer &= mask;
        cur_inv_mmer >>= 2;
        cur_inv_mmer |= ( (0x2 ^ encoded) << (2 * (m - 1))); 
        
        if (canonical_m_mers){
            canon  = (current_mmer < cur_inv_mmer) ? current_mmer : cur_inv_mmer;
        } else {
            canon  = current_mmer;
        }

        const uint64_t s_hash = MinimizerHasher::hash(canon, seed);
        minv                  = (s_hash < minv) ? s_hash : minv;

        if( minv == buffer[0] ){
            minv = s_hash;
            for(int p = 0; p < z; p += 1) {
                const uint64_t value = buffer[p + 1];
                minv = (minv < value) ? minv : value;
                buffer[p] = value;
            }
            buffer[z] = s_hash; // on memorise le hash du dernier m-mer
        }else{
            for(int p = 0; p < z; p += 1) {
                buffer[p] = buffer[p+1]; //TODO: pas moyen d'eviter les deplacements de donnée ici ?
            }
            buffer[z] = s_hash; // on memorise le hash du dernier m-mer
        }


        if( accumulator[n_minmer-1].hash != minv ){
            accumulator.push_back({minv, 1});
            n_minmer += 1;
        } else {
            accumulator[n_minmer-1].size += 1;
        }

        cnt += 1;
    }

    mm_count = n_minmer;
    return contig_size - k + 1;
}

}  // namespace minimizer

#endif // MINIMIZER_HPP