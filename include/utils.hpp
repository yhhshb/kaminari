#ifndef KAMINARI_UTILS_HPP
#define KAMINARI_UTILS_HPP

#include <vector>
#include <string>
#include <cstdint>
#include <iostream>

#include "sys/types.h"
#include "sys/sysinfo.h"

#include "stdlib.h"
#include "stdio.h"
#include "string.h"

namespace kaminari::utils {

std::vector<std::string> read_filenames(std::string const& filenames_list);
std::string get_tmp_filename(const std::string& tmp_dirname, const std::string& prefix, uint64_t run_identifier);

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

int parseLine(char* line);
uint64_t getTotalVirtualMem();
uint64_t getVirtualMemUsed();
int getVirtualMemUsedByProcess();
uint64_t getTotalPhysMem();
uint64_t getPhysMemUsed();
int getPhysMemUsedByProcess();
void printRAMInfo();

} // namespace kaminari::utils

#endif // KAMINARI_UTILS_HPP