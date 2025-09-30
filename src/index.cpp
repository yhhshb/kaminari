#include "../include/index.hpp"

namespace kaminari {
namespace minimizer {

index::index()
    :
    nb_docs(0),
    k(0),
    m(0),
    b(0),
    seed(0),
    canonical(false),
    pthash_constant(0)
{}

index::index(build::options_t& build_parameters)
    : 
    nb_docs(build_parameters.input_filenames.size()),
    k(build_parameters.k),
    m(build_parameters.m),
    b(build_parameters.b),
    seed(build_parameters.seed),
    canonical(build_parameters.canonical),
    pthash_constant(build_parameters.pthash_constant)
{
    build(build_parameters); //see index_build.hpp
}


void index::memory_breakdown(std::ostream& out) const noexcept
{
    libra scale;
    out << "The MPHF of minimizers weights: " << hf.num_bits() / 8 << " Bytes\n";
    //scale.visit(m_ccs);
    //out << "Colors weight: " << scale.get_byte_size() << " Bytes\n";
    //scale.reset();
    //TODO:fix this
    //scale.visit(m_map);
    //out << "The mapping from minimizers to colors weights: " << m_map.num_bytes() << " Bytes\n";
    //scale.reset();
}


typename index::pthash_opt_t index::get_pthash_options(build::options_t& build_parameters)
{
    pthash_opt_t opts;
    opts.seed = build_parameters.seed;
    opts.lambda = build_parameters.pthash_constant; // (too slow = try decreasing), higher lambda : more space efficient 
    opts.alpha = 0.97; //was 0.94
    opts.search = pthash::pthash_search_type::add_displacement;
    opts.avg_partition_size = 3000;
    opts.verbose = (build_parameters.verbose > 0);
    
    opts.ram = build_parameters.max_ram_MB * constants::GB;
    opts.num_threads = build_parameters.nthreads;
    opts.tmp_dir = build_parameters.output_dirname + "/tmp";

    opts.dense_partitioning = true;

    return opts;
}

} // namespace minimizer
} // namespace kaminari