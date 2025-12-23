#include <string>
#include <iostream>
#include <fstream>

#include "../include/constants.hpp"
#include "../include/index.hpp"
#include "../include/build/index_build.hpp"
#include "../include/build/build_options.hpp"
#include "../include/build/build.hpp"
#include "../include/compact_vector.hpp"
#include "../include/utils.hpp"

namespace kaminari::build {

options_t check_args(const argparse::ArgumentParser& parser);

int main(const argparse::ArgumentParser& parser) 
{
    auto start_time = std::chrono::high_resolution_clock::now();

    auto opts = check_args(parser);
    utils::create_directory(opts.output_dirname);
    utils::create_directory(opts.output_dirname+"/tmp");
    minimizer::index idx(opts); //builds colormapper (CM) & colorsets (CS)

    std::ofstream out(opts.output_dirname + "/index.kaminari", std::ios::binary);
    saver saver(out);
    idx.visit(saver); //only MPHF & metadata; CM & CS saved during index build

    if (opts.verbose >= 2) std::cerr << "[II] Written " << saver.get_byte_size() << " Bytes of metadata\n";

    if (opts.verbose >= 1){
        std::cout << "[I] Build time: " <<
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - start_time).count() <<
            " milliseconds\n";
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
    parser.add_argument("-o", "--output-dirname")
        .help("output index directory name")
        .default_value("kaminari_index");
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
    parser.add_argument("-b", "--bit-check")
        .help("number of bits used to check minmers")
        .scan<'d', size_t>()
        .default_value(size_t(1));
    parser.add_argument("-g", "--max-ram")
        .help("RAM limit (MB)")
        .scan<'d', std::size_t>()
        .default_value(std::size_t(8192));
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
    parser.add_argument("-K", "--keep-tmp-files")
        .help("keep temporary files")
        .default_value(false)
        .implicit_value(true);
    parser.add_argument("-v", "--verbose")
        .help("tool verbosity level (0 = silent, 1 = times, 2 = times & stats, 3 = debug)")
        .scan<'d', std::size_t>()
        .default_value(std::size_t(0));
    return parser;
}

options_t check_args(const argparse::ArgumentParser& parser)
{
    options_t opts;
    std::size_t tmp;
    opts.input_filenames = parser.get<std::vector<std::string>>("-i");
    opts.output_dirname = parser.get<std::string>("-o");
    
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
        std::cerr << "Warning: max RAM = 0, setting it to 1 GB\n";
        tmp = 1024;
    }
    opts.max_ram_MB = tmp;
    opts.seed = parser.get<uint64_t>("--seed");
    opts.b = parser.get<std::size_t>("--bit-check");
    opts.pthash_constant = parser.get<double>("--pthash-constant");
    opts.canonical = parser.get<bool>("--canonical");
    opts.keep_tmp_files = parser.get<bool>("--keep-tmp-files");
    opts.verbose = parser.get<std::size_t>("--verbose");

    if (opts.input_filenames.size() == 1){
        opts.input_filenames = utils::read_filenames(opts.input_filenames.at(0));
    }
    return opts;
}

} // namespace kaminari