#include <string>
#include <iostream>
#include <fstream>
#include <thread>

#include <chrono>
#include "../include/query.hpp"


namespace kaminari::query {

options_t check_args(const argparse::ArgumentParser& parser);

int main(const argparse::ArgumentParser& parser) 
{
    auto opts = check_args(parser);

    //load index in mem
    minimizer::index<color_classes::hybrid, pthash::compact_vector> idx;
    {
        std::ifstream in(opts.index_filename, std::ios::binary);
        loader loader(in);
        idx.visit(loader);
        if (opts.verbose) std::cerr << "Read " << loader.get_byte_size() << " Bytes\n";
    }

    //open output file
    std::streambuf * buf;
    std::ofstream of;
    if(opts.output_filename != "") {
        of.open(opts.output_filename.c_str());
        buf = of.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }
    std::ostream out(buf);


    std::chrono::steady_clock::time_point begin;
    begin = std::chrono::steady_clock::now();

    if (opts.nthreads == 1) {
        opts.nthreads += 1;
        std::cerr << "1 thread was specified, but an additional thread will be allocated for parsing\n";
    }
    
    //prep parser + threads
    fastx_parser::FastxParser<fastx_parser::ReadSeq> rparser(opts.input_filenames, opts.nthreads, opts.nthreads - 1);

    rparser.start();
    std::vector<std::thread> workers;
    std::mutex ofile_mut;

    if (opts.ranking){
        if (opts.verbose > 0){
            std::cerr << "Start queries, ranking=on, threshold=" << opts.threshold_ratio << "\n";
        } 
        for (uint64_t i = 1; i != opts.nthreads; ++i) {
            workers.push_back(
                std::thread(
                    [&idx, &rparser, &opts, &out, &ofile_mut]() {
                        ranking_queries(idx, rparser, opts, out, ofile_mut);
                    }
                )
            );
        }
    }
    else{
        if (opts.verbose > 0){
            std::cerr << "Start queries, ranking=off, threshold=" << opts.threshold_ratio << "\n";
        } 
        for (uint64_t i = 1; i != opts.nthreads; ++i) {
            workers.push_back(
                std::thread(
                    [&idx, &rparser, &opts, &out, &ofile_mut]() {
                        classic_queries(idx, rparser, opts, out, ofile_mut);
                    }
                )
            );
        }
    }
    
    for (auto& w : workers) { w.join(); }
    rparser.stop();

    std::cerr << "Query all seqs in " << opts.input_filenames << " took " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - begin).count() << "ms" << std::endl;

    return 0;
}

void ranking_queries(minimizer::index<color_classes::hybrid, pthash::compact_vector>& idx, fastx_parser::FastxParser<fastx_parser::ReadSeq>& rparser, options_t& opts, std::ostream& outstream, std::mutex& ofile_mut){
    auto rg = rparser.getReadGroup();
    std::stringstream ss;
    uint64_t buff_size = 0;
    constexpr uint64_t buff_thresh = 50; //Hardcoded buffer size
    
    while (rparser.refill(rg)) {
        for (auto const& record : rg) {
            auto ids = idx.ranking_query_union_threshold(record.seq.c_str(), record.seq.size(), opts);
            std::sort(
                ids.begin(),
                ids.end(),
                [](auto const& x, auto const& y) { return x.score > y.score; }
            );
            buff_size += 1;

            ss << record.name << "\t" << ids.size();
            for (auto c : ids) { ss << "\t(" << c.item << "," << c.score << ")"; }
            ss << "\n";

            if (buff_size > buff_thresh) {
                std::string outs = ss.str();
                ss.str("");
                ofile_mut.lock();
                outstream.write(outs.data(), outs.size());
                ofile_mut.unlock();
                buff_size = 0;
            }
        }
    }
    // dump anything left in the buffer
    if (buff_size > 0) {
        std::string outs = ss.str();
        ss.str("");
        ofile_mut.lock();
        outstream.write(outs.data(), outs.size());
        ofile_mut.unlock();
        buff_size = 0;
    }
}


void classic_queries(minimizer::index<color_classes::hybrid, pthash::compact_vector>& idx, fastx_parser::FastxParser<fastx_parser::ReadSeq>& rparser, options_t& opts, std::ostream& outstream, std::mutex& ofile_mut){
    using query_fn_t = std::vector<uint32_t> (minimizer::index<color_classes::hybrid, pthash::compact_vector>::*)(const char*, std::size_t, options_t&) const;

   query_fn_t query_algo;
    if (opts.threshold_ratio == 1) { //TODO: remove this, it is a temporary fix (full intersec deleted)
        query_algo = &minimizer::index<color_classes::hybrid, pthash::compact_vector>::query_union_threshold;
    } else {
        query_algo = &minimizer::index<color_classes::hybrid, pthash::compact_vector>::query_union_threshold;
    }

   auto rg = rparser.getReadGroup();
   std::stringstream ss;
   uint64_t buff_size = 0;
   constexpr uint64_t buff_thresh = 50; //Hardcoded buffer size

   while (rparser.refill(rg)) {
        for (auto const& record : rg) {
            auto ids = (idx.*query_algo)(record.seq.c_str(), record.seq.size(), opts);
            buff_size += 1;

            ss << record.name << "\t" << ids.size();
            for (auto c : ids) { ss << "\t" << c; }
            ss << "\n";

            if (buff_size > buff_thresh) {
                std::string outs = ss.str();
                ss.str("");
                ofile_mut.lock();
                outstream.write(outs.data(), outs.size());
                ofile_mut.unlock();
                buff_size = 0;
            }
        }
    }
    // dump anything left in the buffer
    if (buff_size > 0) {
        std::string outs = ss.str();
        ss.str("");
        ofile_mut.lock();
        outstream.write(outs.data(), outs.size());
        ofile_mut.unlock();
        buff_size = 0;
    }
}


argparse::ArgumentParser get_parser()
{
    argparse::ArgumentParser parser("query");
    parser.add_description("Query a kaminari index");
    parser.add_argument("-x", "--index")
        .help("kaminari index to use")
        .required();
    parser.add_argument("-i", "--input-list")
        .help("list of fasta files to query, if 1 file ending with \".list\" is provided, it is assumed to be a list of filenames")
        .nargs(argparse::nargs_pattern::at_least_one)
        .required();
    parser.add_argument("-o", "--output-filename")
        .help("query output, for each sequence of each input file, output its color")
        .default_value("");
    parser.add_argument("-d", "--tmp-dir")
        .help("temporary directory")
        .default_value(std::string("."));
    parser.add_argument("-r", "--ratio")
        .help("ratio of kmer needed to select a color (e.g. r=0.3 -> need atleast 30\% of kmers belonging to the color c1 to select c1)")
        .default_value(std::string("1.0"));
    parser.add_argument("--ranking") 
        .help("flag - list of document ids will be ranked according to the number of kmers/doc [default: off]")
        .default_value(false)   // Default is off
        .implicit_value(true);   // Enabled when provided
    parser.add_argument("-g", "--max-ram")
        .help("RAM limit (GB) [4]")
        .scan<'d', std::size_t>()
        .default_value(std::size_t(4));
    parser.add_argument("-t", "--threads")
        .help("number of threads")
        .scan<'u', std::size_t>()
        .default_value(std::size_t(1));
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

    opts.ranking = parser.get<bool>("--ranking");
    

    opts.verbose = parser.get<std::size_t>("--verbose");

    // if (opts.check and opts.nthreads != 1) {
    //     std::cerr << "[Warning] Checking does not support multi-threading\n";
    // }

    // check if input is a list of filenames or fasta files
    if (opts.input_filenames.size() == 1 && 
        opts.input_filenames.at(0).size() >= 5 &&
        opts.input_filenames.at(0).substr(opts.input_filenames.at(0).size() - 5) == ".list") 
    {
        opts.input_filenames = utils::read_filenames(opts.input_filenames.at(0));
    }

    return opts;
}

} // namespace kaminari