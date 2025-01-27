#include <iostream>
#include <cinttypes>
#include <random>

#include "bundled/biolib/include/hash.hpp"
#include "bundled/pthash/include/pthash.hpp"

#define P10_UINT64 10000000000000000000ULL   /* 19 zeroes */
#define E10_UINT64 19

#define STRINGIZER(x)   # x
#define TO_STRING(x)    STRINGIZER(x)

typedef unsigned __int128 uint128_t;

typedef pthash::build_configuration pthash_opt_t;
typedef pthash::dense_partitioned_phf<pthash::murmurhash2_64,               pthash::opt_bucketer, pthash::mono_EF, true, pthash::pthash_search_type::add_displacement> pthash_minimizers_mphf_t;

static int print_u128_u(uint128_t u128)
{
    int rc;
    if (u128 > UINT64_MAX)
    {
        uint128_t leading  = u128 / P10_UINT64;
        uint64_t  trailing = u128 % P10_UINT64;
        rc = print_u128_u(leading);
        rc += printf("%." TO_STRING(E10_UINT64) PRIu64, trailing);
    }
    else
    {
        uint64_t u64 = u128;
        rc = printf("%" PRIu64, u64);
    }
    return rc;
}


pthash_opt_t get_pthash_options(){
    pthash_opt_t opts;
    opts.seed = 42;
    opts.lambda = 4; // (too slow = try decreasing), higher lambda : more space efficient 
    opts.alpha = 0.97; //was 0.94
    opts.search = pthash::pthash_search_type::add_displacement;
    opts.avg_partition_size = 3000;
    opts.verbose = 2;
    
    opts.ram = 256ULL * 1000000000;
    opts.num_threads = 32;
    opts.tmp_dir = "/tmp";

    opts.dense_partitioning = true;

    return opts;
}


int main() {
    hash::double_hash64 hash;
    auto res = hash("hello", 5);

    __uint128_t hash_128 = (static_cast<__uint128_t>(res[1]) << 64) | res[0];
    //print_u128_u(hash_128); 


    std::vector<__uint128_t> minmers;

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dis;
    for (size_t i = 0; i < 2600000000; ++i) {
        __uint128_t high = dis(gen);
        __uint128_t low = dis(gen);
        minmers.push_back((high << 64) | low);
    }



    pthash_minimizers_mphf_t hf;
    hf.build_in_internal_memory(minmers.begin(), minmers.size(), get_pthash_options());

    return 0;
}