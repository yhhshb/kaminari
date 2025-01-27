#include <cstdint>
#include <sstream>
#include <fstream>
#include "../include/utils.hpp"

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

// std::string get_tmp_filename(const std::string& prefix, uint16_t batch_id, uint16_t depth, std::thread::id tid) // FIXME thread::id is not unique, it can be reused by new threads
// {
//     std::ostringstream oss;
//     oss << tid;  // Convert thread ID to a string representation
//     std::string thread_id_str = oss.str();

//     std::stringstream filename;
//     filename << prefix << "_" << batch_id << "_" << depth << "_" << thread_id_str << ".bin";
//     return filename.str();
// }

} // namespace util