#pragma once

#include <argparse/argparse.hpp>
#include "query_options.hpp"
#include "index.hpp"
#include "hybrid.hpp"
#include "utils.hpp"
#include "constants.hpp"
#include "query_options.hpp"
#include "../bundled/biolib/bundled/prettyprint.hpp"

namespace kaminari::query {

void ranking_queries(minimizer::index<color_classes::hybrid, pthash::compact_vector>& idx, kseq_t* seq, std::ostream& out, options_t& opts);

void classic_queries(minimizer::index<color_classes::hybrid, pthash::compact_vector>& idx, kseq_t* seq, std::ostream& out, options_t& opts);

argparse::ArgumentParser get_parser();
int main(const argparse::ArgumentParser& parser);

} // namespace kaminari