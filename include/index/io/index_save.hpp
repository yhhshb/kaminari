#ifndef INDEX_SAVE_HPP
#define INDEX_SAVE_HPP

#include <filesystem>

namespace kaminari {
namespace minimizer {

#define CLASS_HEADER template <class ColorClasses, class ColorMapper>
#define METHOD_HEADER index<ColorClasses, ColorMapper>



CLASS_HEADER
size_t METHOD_HEADER::save(std::string const& dirname)
{
    std::printf("Saving index to disk...\n");
    std::size_t bytes_written = 0;

    // Remove the directory if it exists
    if (std::filesystem::exists(dirname)) {
        std::filesystem::remove_all(dirname);
    }

    // Create the index directory
    std::filesystem::create_directories(dirname);

    //save core (everything except ColorMapper)
    std::string index_file = dirname + "/index.kaminari";
    std::ofstream out(index_file, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Failed to open output file: " + index_file);
    }
    saver saver(out);
    visit(saver);

    // Create a subdirectory for the mapper
    std::filesystem::create_directories(dirname + "/mapper");


    bytes_written = saver.get_byte_size();
    return bytes_written;
}


#undef CLASS_HEADER
#undef METHOD_HEADER


} // namespace minimizer
} // namespace kaminari


#endif // INDEX_BUILD_HPP