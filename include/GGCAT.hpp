#ifndef GGCAT_HPP
#define GGCAT_HPP

#include <functional>
#include <vector>

#include "../bundled/ggcat/crates/capi/ggcat-cpp-api/include/ggcat.hh"
#include "constants.hpp"

namespace kaminari {

class GGCAT 
{
    public:
        GGCAT(opt_t const& options);
        void loop_through_unitigs(std::function<void(ggcat::Slice<char> const, ggcat::Slice<uint32_t> const, bool)> callback, uint64_t num_threads = 1) const;
        ~GGCAT();

    private:
        void build(opt_t::fn_t const& filenames_list, uint64_t mem_gigas, uint64_t k, uint64_t num_threads, std::string const& tmp_dirname, std::string const& output_basename);
        ggcat::GGCATInstance* m_instance;
        std::string m_graph_file;
        std::string m_output_filename;
        uint64_t m_k;
};

}  // namespace fulgor

#endif // GGCAT_HPP