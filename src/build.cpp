#include <string>
#include <iostream>
#include <fstream>

#include "../include/constants.hpp"
#include "../include/index/index.hpp"
#include "../include/index/index_build.hpp"
#include "../include/index/index_build_metage.hpp"
#include "../include/hybrid.hpp"
#include "../include/utils.hpp"
#include "../include/build_options.hpp"
#include "../include/build.hpp"
#include "../include/compact_vector.hpp"

namespace kaminari::build {

options_t check_args(const argparse::ArgumentParser& parser);

int main(const argparse::ArgumentParser& parser) 
{
    auto opts = check_args(parser);
    minimizer::index<color_classes::hybrid, kaminari::compact_vector> idx(opts);
    if (opts.verbose) {
        idx.memory_breakdown(std::cerr);
        std::cerr << "\n";
    }
    if (opts.output_filename != "") {
        std::ofstream out(opts.output_filename, std::ios::binary);
        saver saver(out);
        idx.visit(saver);
        if (opts.verbose) std::cerr << "Written " << saver.get_byte_size() << " Bytes\n";
    }
    return 0;
}

argparse::ArgumentParser get_parser()
{
    argparse::ArgumentParser parser("build", "1.0.0", argparse::default_arguments::help);
    parser.add_description("Build a kaminari index from a list of datasets (build -h for more information)");

    parser.add_argument("-i", "--input-list")
        .help("list of input files")
        .nargs(argparse::nargs_pattern::at_least_one)
        .required();
    parser.add_argument("-o", "--output-filename")
        .help("output index filename (\".kaminari\" advised)")
        .default_value("index.kaminari");
    parser.add_argument("-k")
        .help("k-mer length")
        .scan<'u', std::size_t>()
        .default_value(size_t(31));
    parser.add_argument("-m")
        .help("minimizer length (must be < k)")
        .scan<'u', std::size_t>()
        .default_value(size_t(19));
    parser.add_argument("-a", "--canonical")
        .help("canonical minimizers")
        .default_value(false)
        .implicit_value(true);
    parser.add_argument("--metagenome")
        .help("recommended for big (>1000docs) metagenomic data collection")
        .default_value(false)
        .implicit_value(true);
    parser.add_argument("-b", "--bit-check")
        .help("number of bits used to check minmers")
        .scan<'d', size_t>()
        .default_value(size_t(1));
    parser.add_argument("-d", "--tmp-dir")
        .help("temporary directory")
        .default_value(std::string("."));
    parser.add_argument("-g", "--max-ram")
        .help("RAM limit (GB)")
        .scan<'d', std::size_t>()
        .default_value(std::size_t(4));
    parser.add_argument("-t", "--threads")
        .help("number of threads")
        .scan<'u', std::size_t>()
        .default_value(std::size_t(1));
    parser.add_argument("-s", "--seed")
        .help("murmurhash random seed")
        .scan<'d', uint64_t>()
        .default_value(uint64_t(42));
    parser.add_argument("-c", "--pthash-constant")
        .help("PTHash build constant, higher = slower but more space efficient")
        .scan<'f', double>()
        .default_value(double(4));
    parser.add_argument("-v", "--verbose")
        .help("increase output verbosity")
        .scan<'d', std::size_t>()
        .default_value(std::size_t(0));
    return parser;
}

options_t check_args(const argparse::ArgumentParser& parser)
{
    options_t opts;
    std::size_t tmp;
    opts.input_filenames = parser.get<std::vector<std::string>>("-i");
    opts.output_filename = parser.get<std::string>("-o");
    opts.tmp_dir = parser.get<std::string>("--tmp-dir");
    
    tmp = parser.get<std::size_t>("-k");
    if (tmp >= constants::MAX_KMER_SIZE) throw std::invalid_argument("k-mer size must be < " + std::to_string(constants::MAX_KMER_SIZE));
    opts.k = static_cast<decltype(opts.k)>(tmp);

    tmp = parser.get<std::size_t>("-m");
    if (tmp > opts.k) throw std::invalid_argument("minimizer size must be < k");
    opts.m = static_cast<decltype(opts.m)>(tmp);

    tmp = std::min(parser.get<std::size_t>("-t"), std::size_t(std::thread::hardware_concurrency()));
    opts.nthreads = static_cast<decltype(opts.nthreads)> (tmp);

    tmp = parser.get<std::size_t>("--max-ram");
    if (tmp == 0) {
        std::cerr << "Warning: max ram = 0, setting it to 1 GB\n";
        tmp = 1;
    }
    opts.max_ram = tmp;
    opts.seed = parser.get<uint64_t>("--seed");
    opts.b = parser.get<std::size_t>("--bit-check");
    opts.pthash_constant = parser.get<double>("--pthash-constant");
    opts.canonical = parser.get<bool>("--canonical");
    opts.metagenome = parser.get<bool>("--metagenome");
    opts.verbose = parser.get<std::size_t>("--verbose");

    if (opts.input_filenames.size() == 1){
        opts.input_filenames = utils::read_filenames(opts.input_filenames.at(0));
    }
    return opts;
}

} // namespace kaminari