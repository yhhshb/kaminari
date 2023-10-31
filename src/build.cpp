#include <string>
#include <iostream>
#include <fstream>
#include <thread>
#include "../include/constants.hpp"
#include "../include/build.hpp"
#include "../include/index.hpp"
#include "../include/hybrid.hpp"
#include "../include/array_mapper.hpp"

namespace kaminari {

opt_t check_args(const argparse::ArgumentParser& parser);

int build_main(const argparse::ArgumentParser& parser) 
{
    auto opts = check_args(parser);
    index<color_classes::hybrid, mapper::array_based<lphash_mphf_t>> idx(opts);
    if (opts.verbose) {
        // idx.print_map_histogram(std::cerr);
        idx.memory_breakdown(std::cerr);
        std::cerr << "\n";
    }
    if (opts.output_filename != "") {
        std::ofstream out(opts.output_filename, std::ios::binary);
        saver saver(out);
        idx.visit(saver);
        if (opts.verbose) std::cerr << "Written " + std::to_string(saver.get_byte_size()) + " Bytes\n";
    }
    return 0;
}

argparse::ArgumentParser get_parser_build()
{
    argparse::ArgumentParser parser("build");
    parser.add_description("Build a kaminari index from a colored compacted de Bruijn Graph");
    parser.add_argument("-i", "--input-list")
        .help("list of input files")
        .nargs(argparse::nargs_pattern::at_least_one)
        .required();
    parser.add_argument("-o", "--output-filename")
        .help("output Index")
        .default_value("");
    parser.add_argument("-k")
        .help("k-mer length")
        .scan<'u', std::size_t>()
        .required();
    parser.add_argument("-m")
        .help("minimizer length (must be < k)")
        .scan<'u', std::size_t>()
        .required();
    parser.add_argument("-t", "--threads")
        .help("number of threads")
        .scan<'u', std::size_t>()
        .default_value(std::size_t(1));
    parser.add_argument("-d", "--tmp-dir")
        .help("temporary directory")
        .default_value(std::string("."));
    parser.add_argument("-g", "--max-ram")
        .help("RAM limit (GB) [4]")
        .scan<'d', std::size_t>()
        .default_value(std::size_t(4));
    parser.add_argument("-C", "--check")
        .help("check MPHF correctness")
        .implicit_value(true)
        .default_value(false);
    parser.add_argument("-v", "--verbose")
        .help("increase output verbosity")
        .default_value(false)
        .implicit_value(true);
    return parser;
}

opt_t::fn_t read_filenames(std::string const& filenames_list) 
{
    opt_t::fn_t buffer;
    std::ifstream in(filenames_list);
    if (!in.is_open()) throw std::runtime_error("error in opening file");
    std::string filename;
    while (in >> filename) buffer.push_back(filename);
    in.close();
    return buffer;
}

opt_t check_args(const argparse::ArgumentParser& parser)
{
    opt_t opts;
    std::size_t tmp;
    opts.input_filenames = parser.get<opt_t::fn_t>("-i");
    opts.output_filename = parser.get<std::string>("-o");
    opts.tmp_dir = parser.get<std::string>("--tmp-dir");
    
    tmp = parser.get<std::size_t>("-k");
    if (tmp >= constants::MAX_KMER_SIZE) throw std::invalid_argument("k-mer size must be < " + std::to_string(constants::MAX_KMER_SIZE));
    opts.k = static_cast<decltype(opts.k)>(tmp);

    tmp = parser.get<std::size_t>("-m");
    if (tmp > opts.k) throw std::invalid_argument("minimizer size must be < k");
    opts.m = static_cast<decltype(opts.m)>(tmp);

    tmp = std::min(parser.get<std::size_t>("-t"), std::size_t(std::thread::hardware_concurrency()));
    opts.nthreads = static_cast<decltype(opts.nthreads) > (tmp);

    tmp = parser.get<std::size_t>("--max-ram");
    if (tmp == 0) {
        std::cerr << "Warning: max ram = 0, setting it to 1 GB\n";
        tmp = 1;
    }
    opts.max_ram = tmp;

    opts.check = parser.get<bool>("--check");
    opts.verbose = parser.get<bool>("--verbose");

    if (opts.check and opts.nthreads != 1) {
        std::cerr << "[Warning] Checking does not support multi-threading\n";
    }
    if (opts.input_filenames.size() == 1) opts.input_filenames = read_filenames(opts.input_filenames.at(0));

    return opts;
}

} // namespace kaminari