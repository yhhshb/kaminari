#include <string>
#include <iostream>
#include <fstream>
#include <thread>
#include "../include/index.hpp"
#include "../include/hybrid.hpp"
#include "../include/utils.hpp"
#include "../include/constants.hpp"
#include "../include/query_options.hpp"
#include "../include/query.hpp"
#include "../bundled/biolib/bundled/prettyprint.hpp"

namespace kaminari::query {

options_t check_args(const argparse::ArgumentParser& parser);

int main(const argparse::ArgumentParser& parser) 
{
    auto opts = check_args(parser);
    minimizer::index<color_classes::hybrid, pthash::compact_vector> idx;
    {
        std::ifstream in(opts.index_filename, std::ios::binary);
        loader loader(in);
        idx.visit(loader);
        if (opts.verbose) std::cerr << "Read " << loader.get_byte_size() << " Bytes\n";
    }

    std::streambuf * buf;
    std::ofstream of;
    if(opts.output_filename != "") {
        of.open(opts.output_filename.c_str());
        buf = of.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }
    std::ostream out(buf);

    gzFile fp = nullptr;
    kseq_t* seq = nullptr;
    for (auto filename : opts.input_filenames) {
        if ((fp = gzopen(filename.c_str(), "r")) == NULL)
            throw std::runtime_error("Unable to open input file " + filename);
        seq = kseq_init(fp);
        while (kseq_read(seq) >= 0) {
            auto ids = idx.query_full_intersection(seq->seq.s, seq->seq.l, opts.threshold_ratio, opts.verbose);
            out << std::string(seq->name.s, seq->name.l) << "\n";
            out << ids << "\n";
        }
        seq = nullptr;
        fp = nullptr;
    }
    return 0;
}

argparse::ArgumentParser get_parser()
{
    argparse::ArgumentParser parser("query");
    parser.add_description("Query a kaminari index");
    parser.add_argument("-x", "--index")
        .help("kaminari index to use")
        .required();
    parser.add_argument("-i", "--input-list")
        .help("list of input files to query")
        .nargs(argparse::nargs_pattern::at_least_one)
        .required();
    parser.add_argument("-o", "--output-filename")
        .help("query output, for each sequence of each input file, output its color")
        .default_value("");
    parser.add_argument("-t", "--threads")
        .help("number of threads")
        .scan<'u', std::size_t>()
        .default_value(std::size_t(1));
    parser.add_argument("-d", "--tmp-dir")
        .help("temporary directory")
        .default_value(std::string("."));
    parser.add_argument("-r", "--ratio")
        .help("ratio of kmer needed to select a color (e.g. r=0.3 -> need atleast 30\% of kmers belonging to the color c1 to select c1)")
        .default_value(std::string("1.0"));
    parser.add_argument("-g", "--max-ram")
        .help("RAM limit (GB) [4]")
        .scan<'d', std::size_t>()
        .default_value(std::size_t(4));
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
    opts.index_filename = parser.get<std::string>("-x");
    opts.input_filenames = parser.get<std::vector<std::string>>("-i");
    opts.output_filename = parser.get<std::string>("-o");
    opts.tmp_dir = parser.get<std::string>("--tmp-dir");

    tmp = std::min(parser.get<std::size_t>("-t"), std::size_t(std::thread::hardware_concurrency()));
    opts.nthreads = static_cast<decltype(opts.nthreads)> (tmp);

    tmp = parser.get<std::size_t>("--max-ram");
    if (tmp == 0) {
        std::cerr << "Warning: max ram = 0, setting it to 1 GB\n";
        tmp = 1;
    }
    opts.max_ram = tmp;

    std::string ratio_str = parser.get<std::string>("--ratio");
    opts.threshold_ratio = std::stof(ratio_str);

    if (opts.threshold_ratio <= 0.0 || opts.threshold_ratio > 1.0) {
        std::cerr << "Warning: kmer ratio needs to be somewhere in ]0.0, 1.0] , setting it to 1.0\n";
        opts.threshold_ratio = 1;
    }
    

    opts.verbose = parser.get<std::size_t>("--verbose");

    // if (opts.check and opts.nthreads != 1) {
    //     std::cerr << "[Warning] Checking does not support multi-threading\n";
    // }
    if (opts.input_filenames.size() == 1) opts.input_filenames = utils::read_filenames(opts.input_filenames.at(0));
    return opts;
}

} // namespace kaminari