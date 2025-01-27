#ifndef KAMINARI_COLOR_CLASSES_HYBRID_HPP
#define KAMINARI_COLOR_CLASSES_HYBRID_HPP

#include <vector>
#include "../bundled/biolib/include/bit_vector.hpp"
#include "../bundled/biolib/include/bit_parser.hpp"
#include "../bundled/biolib/include/codes.hpp"
#include "../bundled/biolib/include/packed_vector.hpp"
#include "../bundled/biolib/include/elias_fano.hpp"

namespace kaminari {
namespace color_classes {

class hybrid
{
    public:
        typedef uint32_t color_t;
        typedef bit::vector<uint64_t> bit_vector;
        typedef bit::packed::vector<uint64_t> packed_vector;
        typedef bit::ef::array ef_sequence;
        typedef bit::parser<uint64_t> bit_parser;

        class builder
        {
            public:
                builder(std::size_t number_of_documents, std::size_t verbosity_level = 0);
                void add_color_set(color_t const * const colors, std::size_t list_size);
                void build(hybrid& index);

            private:
                uint32_t m_num_docs;
                uint32_t m_sparse_set_threshold_size;
                uint32_t m_very_dense_set_threshold_size;
                std::size_t m_num_lists;
                std::size_t m_num_total_integers;
                std::size_t m_verbosity_level;

                bit_vector m_bvb;
                std::vector<uint64_t> m_offsets;
        };

        class row_accessor
        {
            public:
                using iterator_category = std::forward_iterator_tag;
                using difference_type   = std::ptrdiff_t;
                using value_type        = color_t;
                using pointer           = value_type*;
                using reference         = value_type&;

                enum list_type { 
                    delta_gaps = 0,
                    bitmap = 1,
                    complementary_delta_gaps = 2
                };

                row_accessor(hybrid const* parent_color_storage, std::size_t start_idx);
                value_type value() const;
                void next();
                void next_geq(uint64_t lower_bound);
                
                void reinit_for_complemented_set_iteration();
                value_type comp_value() const;
                void comp_next();
                void comp_next_geq(uint64_t lower_bound);

                void reset();

                list_type type() const;
                std::size_t size() const;
                
            private:
                hybrid const* m_parent;
                bit_parser m_parser;
                uint64_t m_begin;
                uint64_t m_orig;
                list_type m_type;

                uint32_t m_size;
                uint32_t m_pos_in_list;
                uint32_t m_prev_val;
                uint32_t m_curr_val;
                
                uint32_t m_comp_list_size;
                uint32_t m_pos_in_comp_list;
                uint32_t m_comp_val;
                
                void find_next();
                
        };

        hybrid();
        row_accessor colors_at(std::size_t color_class_id) const;
        std::size_t num_docs() const;
        std::size_t num_color_classes() const;
        std::size_t num_bits() const;
        void print_stats(std::ostream& out) const;

        template <typename Visitor>
        void visit(Visitor& visitor) 
        {
            visitor.visit(m_num_docs);
            visitor.visit(m_sparse_set_threshold_size);
            visitor.visit(m_very_dense_set_threshold_size);
            visitor.visit(m_offsets);
            visitor.visit(m_colors);
        }

        template <typename Visitor>
        void visit(Visitor& visitor) const
        {
            visitor.visit(m_num_docs);
            visitor.visit(m_sparse_set_threshold_size);
            visitor.visit(m_very_dense_set_threshold_size);
            visitor.visit(m_offsets);
            visitor.visit(m_colors);
        }

    private:
        uint32_t m_num_docs;
        uint32_t m_sparse_set_threshold_size;
        uint32_t m_very_dense_set_threshold_size;
        ef_sequence m_offsets;
        bit_vector m_colors;
        // bool verbose;
};

} // namespace color_classes 
} // namespace kaminari

#endif // KAMINARI_COLOR_CLASSES_HYBRID_HPP