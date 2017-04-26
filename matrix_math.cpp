#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <type_traits>
#include <utility>

namespace matrix_math
{
    using index_t = std::size_t;
    using default_T = std::int_fast64_t;

    template <index_t Width, typename T = default_T>
    class row
    {
        using storage_t = std::array<T, Width>;
        storage_t storage{};

        public:
        using type = T;
        using self_type = row<Width, T>;
        
        row() : storage{} { }
        row(const std::initializer_list<T> init) : storage{init} { }
        
        constexpr T& operator[] (const index_t x) { return storage[x]; }
        constexpr const T operator[] (const index_t x) const { return storage[x]; }

        template <typename Callable>
        void apply(Callable f)
        // Apply a function f(element) to every element.
        {
            for (auto& element : storage)
                f(element);
        }

        // Row multiplication by a constant.
        void operator*= (const T rhs) { apply([rhs] (auto& el) {el *= rhs;}); }
        self_type operator* (const T rhs) const
        {
            self_type new_row{};
            for (index_t i = 0; i < Width; ++i)
                new_row[i] = storage[i] * rhs;
            return new_row;
        }

        // Addition of two rows.
        void operator+= (const row<Width,T>& rhs)
        {
            for (index_t i = 0; i < Width; ++i)
                storage[i] += rhs[i];
        }
    };

    template <index_t Height, index_t Width, typename T = default_T>
    class matrix
    {
        public:
        using type = T;
        using row_t = row<Width, T>;

        private:
        using storage_t = std::array<row<Width, T>, Height>;
        storage_t storage{};

        public:
        matrix() : storage{} { }
        matrix(const std::initializer_list<std::initializer_list<T>> init) : storage{init} { }
        
        constexpr row_t& operator[] (const index_t y) { return storage[y]; }
        constexpr const row_t& operator[] (const index_t y) const { return storage[y]; }

        void swap_rows(const index_t a, const index_t b)
        {
            using std::swap;
            swap(storage[a], storage[b]);
        }

        void add_multiply_rows(const index_t a, const index_t b, const T multiple)
        {
            storage[a] += storage[b] * multiple;
        }
    };

    template <index_t Height, index_t WidthLeft, index_t WidthRight, typename T = default_T>
    auto horizontal_concat(const matrix<Height, WidthLeft, T>& lhs, const matrix<Height, WidthRight, T>& rhs)
        -> matrix<Height, WidthLeft+WidthRight, T>
    {
        matrix<Height, WidthLeft+WidthRight, T> concatenation {};
        for (index_t r = 0; r < Height; ++r)
        {
            for (index_t c = 0; c < WidthLeft; ++c)
                concatenation[r][c] = lhs[r][c];
            for (index_t c = 0; c < WidthRight; ++c)
                concatenation[r][c+WidthLeft] = rhs[r][c];
        }
        return concatenation; 
    }

    // Forward declaration.
    template <index_t Size, typename T = default_T>
    class square_matrix;
    
    template <index_t Size, typename T = default_T>
    const square_matrix<Size, T> identity_matrix = []()
    {
        square_matrix<Size, T> identity {};
        for (index_t i = 0; i < Size; ++i)
            identity[i][i] = T(1);
        return identity;
    };

    template <index_t Size, typename T>
    class square_matrix : public matrix<Size, Size, T>
    {
        public:
        void invert()
        {
            auto augmented_matrix = horizontal_concat(this, identity_matrix<Size, T>);
            augmented_matrix.row_reduce();
            swap(this, augmented_matrix.horizontal_slice(Size));
        }
    };
}

using namespace matrix_math;

square_matrix<3> invert(square_matrix<3> mtx) 
{
    mtx.invert();
    return mtx;
}
