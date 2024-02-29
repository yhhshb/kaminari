#ifndef KAMINARI_COLOR_CLASSES_HYBRID_HPP
#define KAMINARI_COLOR_CLASSES_HYBRID_HPP

#include <vector>
#include "constants.hpp"

namespace kaminari {
namespace color_classes {

class hybrid
{
    public:
        enum list_type { 
            delta_gaps = 0,
            bitmap = 1,
            complementary_delta_gaps = 2
        };

        class builder
        {
            public:
                builder(std::size_t number_of_documents, bool verbose = false);
                void add_color_set(uint32_t const * const colors, uint64_t list_size);
                void build(hybrid& index);

            private:
                uint32_t m_num_docs;
                uint32_t m_sparse_set_threshold_size;
                uint32_t m_very_dense_set_threshold_size;
                uint64_t m_num_lists;
                uint64_t m_num_total_integers;
                bool verbose;

                bit_vector m_bvb;
                std::vector<uint64_t> m_offsets;
        };

        class iterator
        {
            public:
                iterator(hybrid const& ptr, uint64_t begin);
                void reinit_for_complemented_set_iteration();
                uint64_t value() const;
                uint64_t comp_value() const;
                uint64_t operator*() const;
                const iterator& operator++();
                void next_comp();
                void next_geq(uint64_t lower_bound);
                std::size_t size() const;
                std::size_t num_docs() const;
                list_type type() const;

            private:
                hybrid const& m_ptr;
                uint64_t m_begin;
                uint32_t m_num_docs;
                list_type m_type;

                bit_parser m_parser;
                uint32_t m_pos_in_list;
                uint32_t m_size;

                uint32_t m_pos_in_comp_list;
                uint32_t m_comp_list_size;

                uint32_t m_comp_val;
                uint32_t m_prev_val;
                uint32_t m_curr_val;
                
                void next_comp_val();
                void next_geq_comp_val(uint64_t lower_bound);
        };

        hybrid();
        iterator colors(uint64_t color_class_id) const;
        std::size_t num_docs() const;
        std::size_t num_color_classes() const;
        std::size_t num_bits() const;
        void print_stats() const;

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
        bool verbose;
};

} // namespace color_classes 
} // namespace kaminari

#endif // KAMINARI_COLOR_CLASSES_HYBRID_HPP