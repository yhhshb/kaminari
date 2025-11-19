#include <sstream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include "include/utils.hpp"

namespace kaminari::utils {

std::vector<std::string>
read_filenames(std::string const& filenames_list) 
{
    std::vector<std::string> buffer;
    std::ifstream in(filenames_list);
    if (!in.is_open()) throw std::runtime_error("error in opening file");
    std::string filename;
    while (in >> filename) buffer.push_back(filename);
    in.close();
    return buffer;
}

std::string get_tmp_filename(const std::string& tmp_dirname, const std::string& prefix, std::size_t run_identifier)
{
    std::stringstream filename;
    filename 
    << tmp_dirname 
    << (tmp_dirname != "" ? "/" : "") 
    << prefix << "_" << run_identifier << ".bin";
    return filename.str();
}

std::string getExecutablePath() 
{
    char result[PATH_MAX];
    ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
    if (count != -1) {
        result[count] = '\0';
        return std::string(result).substr(0, std::string(result).find_last_of("/"));;
    }
    return "";
}

void create_directory(std::string dirname) 
{
    if (dirname.empty()) return;
    struct stat st;
    if (stat(dirname.c_str(), &st) != 0) {
        if (mkdir(dirname.c_str(), 0755) != 0) {
            throw std::runtime_error("Failed to create directory: " + dirname);
        }
    } else if (!S_ISDIR(st.st_mode)) {
        throw std::runtime_error(dirname + " exists and is not a directory");
    }
}

uint64_t bits_needed(uint64_t b) {
    if (b == 0) return 1;
#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
    // Fast hardware instruction (BSR equivalent)
    return 64 - __builtin_clzll(b);
#else
    // Portable fallback
    uint64_t bits = 0;
    do {
        bits++;
        b >>= 1;
    } while (b);
    return bits;
#endif
}

uint64_t sparse_colors_bits(uint64_t total_colors) {
    uint64_t bits = bits_needed(total_colors - 1);    
    if (bits <= 8) return 8;
    else if (bits <= 16) return 16;
    else if (bits <= 32) return 32;
    else return 64;
}


uint64_t get_file_size(const std::string& filen) {
    struct stat file_status;
    if (stat(filen.c_str(), &file_status) < 0) {
        return -1;
    }
    return file_status.st_size;
}

} // namespace utils