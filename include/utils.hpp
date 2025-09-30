#ifndef KAMINARI_UTILS_HPP
#define KAMINARI_UTILS_HPP

#include <vector>
#include <string>
#include <cstdint>
#include <iostream>
#include <thread>
#include <unistd.h>
#include <limits.h>

namespace kaminari::utils {

std::vector<std::string>
read_filenames(std::string const& filenames_list);

std::string get_tmp_filename(const std::string& tmp_dirname, const std::string& prefix, std::size_t run_identifier);

std::string getExecutablePath();

void create_directory(std::string dirname);

//defined here because template
template<class ForwardIt>
ForwardIt unique_accumulate(ForwardIt first, ForwardIt last)
{
    if (first == last)
        return last;
 
    ForwardIt result = first;
    while (++first != last) {
        // if the first elements are the same, accumulate the second values
        if (result->first == first->first) {
            result->second += first->second;
        } else {
            // move to the next unique element
            if (++result != first) {
                *result = std::move(*first);
            }
        }
    }
 
    return ++result;
}


} // namespace util

#endif // KAMINARI_UTILS_HPP