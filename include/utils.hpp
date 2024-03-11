#ifndef KAMINARI_UTILS_HPP
#define KAMINARI_UTILS_HPP

#include <vector>
#include <string>

namespace kaminari::utils {

std::vector<std::string> read_filenames(std::string const& filenames_list);
std::string get_tmp_filename(const std::string& tmp_dirname, const std::string& prefix, uint64_t run_identifier);

} // namespace kaminari::utils

#endif // KAMINARI_UTILS_HPP