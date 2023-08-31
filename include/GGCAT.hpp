#pragma once

#include <functional>
#include <vector>

#include "../bundled/ggcat/crates/capi/ggcat-cpp-api/include/ggcat.hh"

namespace kaminari {

class GGCAT 
{
    public:
        GGCAT();
        ~GGCAT();
        void build(std::string const& filenames_list, uint64_t mem_gigas, uint64_t k, uint64_t num_threads, std::string const& tmp_dirname, std::string const& output_basename);
        void loop_through_unitigs(std::function<void(ggcat::Slice<char> const, ggcat::Slice<uint32_t> const, bool)> callback, uint64_t num_threads = 1) const;
        uint64_t num_docs() const;
        std::vector<std::string> const& filenames() const;

    private:
        ggcat::GGCATInstance* m_instance;
        std::string m_graph_file;
        std::vector<std::string> m_filenames;
        std::string m_output_filename;
        uint64_t m_k;
};

}  // namespace fulgor