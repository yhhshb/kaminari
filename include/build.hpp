#pragma once

#include <argparse/argparse.hpp>

namespace kaminari {

argparse::ArgumentParser get_parser_build();
int build_main(const argparse::ArgumentParser& parser);

} // namespace kaminari