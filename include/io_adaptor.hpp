#pragma once

#include "../bundled/biolib/include/io.hpp"
#include "../bundled/lphash/external/pthash/external/essentials/include/essentials.hpp"

// use biolib's io functionality with pthash's syntax (due to compatibility)

namespace kaminari {

class saver : public essentials::saver {
    public:
        saver(std::ostream& output_stream) : essentials::saver(output_stream) {}

        template <typename T>
        void apply(T& val) {visit(val);}

        size_t bytes() {return get_byte_size();}
};

class loader : public essentials::loader {
    public:
        loader(std::istream& input_stream) : essentials::loader(input_stream) {}

        template <typename T>
        void visit(T& val) {apply(val);}

        size_t bytes() {return get_byte_size();}
};

}