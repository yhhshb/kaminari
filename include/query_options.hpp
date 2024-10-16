#ifndef QUERY_OPTIONS_HPP
#define QUERY_OPTIONS_HPP

#include <vector>
#include <string>
#include <cstdint>

namespace kaminari::query {

struct options_t {
    std::string index_filename;
    std::vector<std::string> input_filenames;
    std::string output_filename;
    std::string tmp_dir;
    uint8_t nthreads;
    std::size_t max_ram;
    float threshold_ratio;
    bool ranking;
    std::size_t verbose;
};

}

#endif