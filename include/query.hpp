#pragma once

#include <argparse/argparse.hpp>
#include "query_options.hpp"
#include "index.hpp"
#include "hybrid.hpp"
#include "utils.hpp"
#include "constants.hpp"
#include "query_options.hpp"
#include "../bundled/biolib/bundled/prettyprint.hpp"
#include "../bundled/FQFeeder/include/FastxParser.hpp"

namespace kaminari::query {

void ranking_queries(minimizer::index<color_classes::hybrid, pthash::compact_vector>& idx, fastx_parser::FastxParser<fastx_parser::ReadSeq>& rparser, options_t& opts, std::ostream& outstream, std::mutex& ofile_mut);

void classic_queries(minimizer::index<color_classes::hybrid, pthash::compact_vector>& idx, fastx_parser::FastxParser<fastx_parser::ReadSeq>& rparser, options_t& opts, std::ostream& out, std::mutex& ofile_mut);

argparse::ArgumentParser get_parser();
int main(const argparse::ArgumentParser& parser);

} // namespace kaminari