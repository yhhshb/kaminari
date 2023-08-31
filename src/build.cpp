#include <string>
#include <iostream>
#include <thread>
#include "../include/constants.hpp"
#include "../include/build.hpp"
#include "../include/index.hpp"

namespace kaminari {

index::opt_t check_args(const argparse::ArgumentParser& parser);

int build_main(const argparse::ArgumentParser& parser) 
{
    auto opts = check_args(parser);
    index idx(opts);
    
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
        .default_value(1);
    parser.add_argument("-d", "--tmp-dir")
        .help("temporary directory")
        .default_value(std::string("."));
    parser.add_argument("-m", "--max-ram")
        .help("RAM limit (MB)")
        .scan<'d', std::size_t>()
        .default_value(std::size_t(1000));
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

index::opt_t check_args(const argparse::ArgumentParser& parser)
{
    index::opt_t opts;
    std::size_t tmp;
    opts.input_filenames = parser.get<index::opt_t::fn_t>("-i");
    opts.output_filename = parser.get<std::string>("-o");
    opts.tmp_dir = parser.get<std::string>("--tmp-dir");
    
    tmp = parser.get<std::size_t>("-k");
    if (tmp >= constants::MAX_KMER_SIZE) throw std::invalid_argument("k-mer size must be < " + std::to_string(constants::MAX_KMER_SIZE));
    opts.k = static_cast<decltype(opts.k)>(tmp);

    tmp = parser.get<std::size_t>("-m");
    if (tmp > opts.k) throw std::invalid_argument("minimizer size must be < k");
    opts.m = static_cast<decltype(opts.m)>(tmp);

    tmp = std::min(parser.get<std::size_t>("-t"), std::size_t(std::thread::hardware_concurrency()));
    opts.nthreads = static_cast<decltype(opts.nthreads)>(tmp);

    tmp = parser.get<std::size_t>("--max-ram");
    opts.max_ram = tmp;

    opts.check = parser.get<bool>("--check");
    opts.verbose = parser.get<bool>("--verbose");

    return opts;
}

} // namespace kaminari