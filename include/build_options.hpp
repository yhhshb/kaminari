#ifndef BUILD_OPTION_HPP
#define BUILD_OPTION_HPP

#include <vector>
#include <string>

namespace kaminari::build {
struct options_t {
    std::vector<std::string> input_filenames;
    uint8_t k;
    uint8_t m;
    uint8_t b;
    uint64_t seed;
    bool canonical;
    double pthash_constant;
    std::string output_filename;
    std::string tmp_dir;
    uint8_t nthreads;
    std::size_t max_ram;
    std::size_t verbose;
};
} // namespace kaminari::build

#endif