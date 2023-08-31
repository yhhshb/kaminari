#ifndef KAMINARI_INDEX_HPP
#define KAMINARI_INDEX_HPP

#include <string>
#include <vector>

namespace kaminari{

class index
{
    public:
        struct opt_t {
            using fn_t = std::vector<std::string>;
            fn_t input_filenames;
            std::string output_filename;
            std::string tmp_dir;
            uint8_t k;
            uint8_t m;
            uint8_t nthreads;
            std::size_t max_ram;
            bool check;
            bool verbose;
        };
        index(const opt_t& build_parameters);

    private:
        // lphash mphf
        // bit vector or not
        // colors
};

} // namespace kaminari

#endif // KAMINARI_INDEX_HPP