#include "../include/GGCAT.hpp"
#include <iostream>
#include <fstream>


namespace kaminari {

GGCAT::GGCAT(opt_t const& options)
{
    build(options.input_filenames, options.max_ram, options.k, options.nthreads, options.tmp_dir, options.output_filename);
}

GGCAT::~GGCAT() 
{
    try {
        std::remove((m_graph_file).c_str()); // ccdbg filename
        std::remove((m_output_filename + ".ggcat.colors.dat").c_str()); // colors tmp filename
    } catch (std::exception const& e) { 
        std::cerr << e.what() << "\n";
    }
}

void 
GGCAT::build(opt_t::fn_t const& filenames_list, uint64_t mem_gigas, uint64_t k, uint64_t num_threads, std::string const& tmp_dirname, std::string const& output_basename)
{
    m_output_filename = output_basename;
    m_graph_file = output_basename + ".ggcat.fa";
    m_k = k;

    ggcat::GGCATConfig config;
    config.use_temp_dir = true;
    config.prefer_memory = true;
    config.intermediate_compression_level = -1;
    config.use_stats_file = false;
    config.stats_file = "";
    config.temp_dir = tmp_dirname;
    config.memory = mem_gigas;
    config.total_threads_count = num_threads;

    // GGCAT bug:
    // This leaks memory (not much) but it can't be easily fixed because
    // this memory is allocated in the Rust API of GGCAT and freed only at the
    // end of the program.
    m_instance = ggcat::GGCATInstance::create(config);

    std::vector<std::string> color_names;
    color_names.reserve(filenames_list.size());
    for (uint64_t i = 0; i != filenames_list.size(); ++i) color_names.push_back(std::to_string(i));

    constexpr bool forward_only = false;
    constexpr bool output_colors = true;
    constexpr size_t min_multiplicity = 1;
    m_instance->build_graph_from_files(
        ggcat::Slice<std::string>(const_cast<std::string*>(filenames_list.data()), // ATTENTION: const_cast
        filenames_list.size()), 
        m_graph_file, 
        m_k, 
        num_threads, 
        forward_only, 
        min_multiplicity, 
        ggcat::ExtraElaborationStep_UnitigLinks, 
        output_colors, 
        ggcat::Slice<std::string>(color_names.data(), color_names.size())
    );
}

void 
GGCAT::loop_through_unitigs(
    std::function<void(
        ggcat::Slice<char> const /* unitig */,
        ggcat::Slice<uint32_t> const /* colors */, 
        bool /* same_color */
    )> callback,
    uint64_t num_threads) const 
{
    if (m_k == 0) throw std::runtime_error("graph must be built first");
    m_instance->dump_unitigs(m_graph_file, m_k, num_threads, num_threads == 1, callback, true);
}

}