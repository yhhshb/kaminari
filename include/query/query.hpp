#pragma once

#include <argparse/argparse.hpp>
#include "query_options.hpp"

namespace kaminari::query {

argparse::ArgumentParser get_parser();
options_t check_args(const argparse::ArgumentParser& parser);
//void ranking_queries(minimizer::index& idx, fastx_parser::FastxParser<fastx_parser::ReadSeq>& rparser, options_t& opts, std::ostream& outstream, std::mutex& ofile_mut);


int main(const argparse::ArgumentParser& parser);

} // namespace kaminari