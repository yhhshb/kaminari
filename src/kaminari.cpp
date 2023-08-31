#include <iostream>
#include "../include/build.hpp"

using namespace kaminari;

int main(int argc, char* argv[])
{
    auto build_parser = get_parser_build();
    argparse::ArgumentParser program(argv[0]);
    program.add_subparser(build_parser);
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << program;
        return 1;
    }
    if (program.is_subcommand_used(build_parser)) return build_main(build_parser);
    else std::cerr << program << std::endl;
    return 0;
}

// debug build: cmake .. -D CMAKE_BUILD_TYPE=Debug -D PIMHASH_USE_SANITIZERS=On