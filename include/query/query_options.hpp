#ifndef QUERY_OPTIONS_HPP
#define QUERY_OPTIONS_HPP

#include <vector>
#include <string>
#include <cstdint>

namespace kaminari::query {

struct options_t {
    std::string index_dirname;
    std::vector<std::string> input_filenames;
    std::string output_filename;
    uint8_t nthreads;
    float threshold_ratio;
    std::size_t verbose;
};

}

#endif