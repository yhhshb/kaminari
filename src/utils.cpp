#include <cstdint>
#include <sstream>
#include <fstream>
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

std::string getExecutablePath() {
    char result[PATH_MAX];
    ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
    if (count != -1) {
        result[count] = '\0';
        return std::string(result).substr(0, std::string(result).find_last_of("/"));;
    }
    return "";
}

} // namespace utils