#include <iostream>
#include <cinttypes>
#include <random>

#include "bundled/biolib/include/hash.hpp"
#include "bundled/biolib/include/external_memory_vector.hpp"
#include "bundled/pthash/include/pthash.hpp"

using namespace pthash;

#define P10_UINT64 10000000000000000000ULL /* 19 zeroes */
#define E10_UINT64 19

#define STRINGIZER(x) #x
#define TO_STRING(x) STRINGIZER(x)

typedef unsigned __int128 uint128_t;

static int print_u128_u(uint128_t u128) {
    int rc;
    if (u128 > UINT64_MAX) {
        uint128_t leading = u128 / P10_UINT64;
        uint64_t trailing = u128 % P10_UINT64;
        rc = print_u128_u(leading);
        rc += printf("%." TO_STRING(E10_UINT64) PRIu64, trailing);
    } else {
        uint64_t u64 = u128;
        rc = printf("%" PRIu64, u64);
    }
    return rc;
}

typedef build_configuration pthash_opt_t;


typedef pthash::partitioned_phf<  
    pthash::xxhash_128, //murmurhash2_64                     
    pthash::opt_bucketer,               
    pthash::compact,                    
    true>
    pthash_minimizers_mphf_t;


pthash_opt_t get_pthash_options() {
    pthash_opt_t opts;
    opts.seed = 42;
    opts.lambda =
        5;  // (too slow = try decreasing), higher lambda : more space efficient
    opts.alpha = 0.94;  // was 0.94
    //opts.search = pthash::pthash_search_type::add_displacement;
    opts.avg_partition_size = 3000000;
    opts.verbose = true;

    opts.ram = 2ULL * 1024*1024*1024;
    opts.num_threads = 1;
    opts.tmp_dir = "/tmp";

    opts.dense_partitioning = true;

    return opts;
}

int main() {
    std::vector<uint64_t> minmers;
    minmers.reserve(800000000);

    for (size_t i = 0; i < 800000000; ++i) { minmers.push_back(i); }

    pthash_minimizers_mphf_t hf;
    hf.build_in_external_memory(minmers.begin(), minmers.size(),
                                get_pthash_options());

    std::cout << "ALL GOOD: taking "
              << static_cast<double>(hf.num_bits()) / hf.num_keys()
              << " bits/key" << std::endl;

    return 0;
}
