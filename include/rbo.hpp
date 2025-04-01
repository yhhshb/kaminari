#pragma once

#include <argparse/argparse.hpp>

namespace kaminari::rbo {

argparse::ArgumentParser get_parser();
int main(const argparse::ArgumentParser& parser);

} // namespace kaminari