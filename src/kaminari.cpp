#include <iostream>
#include "../include/build.hpp"
#include "../include/query.hpp"
#include "../include/rbo.hpp"

using namespace kaminari;

int main(int argc, char* argv[])
{
    argparse::ArgumentParser program("kaminari", "1.0.0", argparse::default_arguments::all);
    auto build_parser = build::get_parser();
    auto query_parser = query::get_parser();
    auto rbo_parser = rbo::get_parser();

    program.add_subparser(build_parser);
    program.add_subparser(query_parser);
    program.add_subparser(rbo_parser);
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << program;
        return 1;
    }
    if (program.is_subcommand_used(build_parser)) return build::main(build_parser);
    else if (program.is_subcommand_used(query_parser)) return query::main(query_parser);
    else if (program.is_subcommand_used(rbo_parser)) return rbo::main(rbo_parser);
    else std::cerr << program << std::endl;
    return 0;
}

// debug build: cmake .. -D CMAKE_BUILD_TYPE=Debug -D PIMHASH_USE_SANITIZERS=On