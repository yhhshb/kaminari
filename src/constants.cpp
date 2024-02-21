#include "../include/constants.hpp"

namespace kaminari::util {

std::string get_tmp_filename(const std::string& tmp_dirname, const std::string& prefix, uint64_t run_identifier)
{
    std::stringstream filename;
    filename << tmp_dirname << "/ii.tmp.run_" << run_identifier << ".bin";
    return filename.str();
}

} // namespace util