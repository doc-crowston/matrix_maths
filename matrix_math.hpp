/*
 * Matrix maths header.
 *
 * Author: Robert H. Crowston, 2017.
 *
 *
 * Requires C++14 or later.
 *
 */

#ifndef CROWSTON_MATRIX_MATH_H
#define CROWSTON_MATRIX_MATH_H

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <ostream>
#include <numeric>
#include <type_traits>
#include <utility>

namespace matrix_math
{
	//
	// Default parameters for matrices.
	// 
    using index_t = std::size_t;
    using default_T = double;
	// This is a fudge factor for floating point types..
	const default_T equality_tolerance{0.00000000001};

	//
	// Exceptions.
	//
	struct matrix_is_degenerate_error : public std::domain_error
	{
		matrix_is_degenerate_error() : std::domain_error("Cannot invert degenerate matrix.") {}
		virtual ~matrix_is_degenerate_error() {}
	};

	//
	// In this system, matrices comprise rows. Access column-by-column is also supported by a
	// special iterator provided by the Matrix class.
	//
	template <index_t Width, typename T = default_T>
    class row
    {
        using storage_t = std::array<T, Width>;
        storage_t storage{};

        public:
        using type = T;
        using self_t = row<Width, T>;
        
		// Constructors.
        row() noexcept : storage{} { }
        row(const std::initializer_list<T> init) noexcept 
		{
			index_t c = 0;
			for (const auto& init_col : init)
				storage[c++] = T(init_col);
		}
        
		// Accessors.
        constexpr T& operator[] (const index_t x) noexcept { return storage[x]; }
        constexpr const T operator[] (const index_t x) const noexcept { return storage[x]; }

		// Iteration (through underlying storage type).
		constexpr auto begin() noexcept { return storage.begin(); }
		constexpr auto end() noexcept { return storage.end(); }

		// Row multiplication by a constant.
        void operator*= (const T rhs) noexcept
		{
			for (auto& element : storage)
				element *= rhs;
		}
        self_t operator* (const T rhs) const noexcept
        {
            self_t new_row{};
            for (index_t i = 0; i < Width; ++i)
                new_row[i] = storage[i] * rhs;
            return new_row;
        }

        // Addition of one row to this one.
        void operator+= (const row<Width,T>& rhs) noexcept
        {
            for (index_t i = 0; i < Width; ++i)
                storage[i] += rhs[i];
        }
    }; // End of class row.

	//
	// Matrix class.
	//
    template <index_t Height, index_t Width, typename T = default_T>
    class matrix
    {
        public:
        using type = T;
        using row_t = row<Width, T>;
		using self_t = matrix<Height, Width, T>;

        private:
        using storage_t = std::array<row<Width, T>, Height>;
        storage_t storage{};

        public:
		// Constructors.
        matrix() noexcept : storage{} { }
        matrix(const std::initializer_list<std::initializer_list<T>> init) noexcept
		{
			index_t r = 0;
			for (const auto& row_init : init)
				storage[r++] = row_t{row_init};
		}

		// Accessors.
        constexpr row_t& operator[] (const index_t y) noexcept { return storage[y]; }
        constexpr const row_t& operator[] (const index_t y) const noexcept 
		{ 
			return storage[y]; 
		}

		// Iteration, by row.
		constexpr auto begin() noexcept { return storage.begin(); }
		constexpr auto end() noexcept   { return storage.end(); }
		
		// Iteration, by column. (Const iterators.)
		constexpr auto column_cbegin(index_t col) const noexcept 
		{ 
			return const_column_iterator(*this, 0, col);
		}
		constexpr auto column_cend(index_t col) const noexcept
		{
			return const_column_iterator(*this, Height, col);
		}

		// Elementary row operations.
        void swap_rows(const index_t a, const index_t b) noexcept
        {
            using std::swap;
            swap(storage[a], storage[b]);
        }

		// The Gauss--Jordan algorithm.
		void row_reduce()
		{
			for (index_t r = 0; r < Height; ++r)
			{
				// We need a non-zero element at position [r][r].
				// Interchange rows to obtain a non-zero.
				if (std::abs(storage[r][r]) <= equality_tolerance)
				{
					index_t s = r+1;
					while (s < Height && std::abs(storage[s][s]) <= equality_tolerance )
						++s;
					if (s == Height)
						throw matrix_is_degenerate_error();
					swap_rows(r,s);
				}
				// We need element at position [r][r] to be 1.
				// Multiply by the reciprocal.
				{
					if (storage[r][r] != T(1))
						storage[r] *= T(1) / storage[r][r];
				}
				// In column r for each row s ≠ r, we need a zero.
				// Subtract from each row s ≠ r, row [r] * [r][s].
				{
					for (index_t s=0; s < Height; ++s)
						if (s != r && storage[s][r] != T(0))
							storage[s] += storage[r] * -storage[s][r];
				}
			}
		}

		// Obtain the right-hand half following row reduction.
		auto get_right_slice() noexcept
			-> matrix<Height, Width/2, T>
		{
			matrix<Height, Width/2, T> right_slice {};
			for (index_t r = 0; r < Height; ++r)
				for (index_t c = 0; c < Width/2; ++c)
					right_slice[r][c] = storage[r][c+Width/2];
			return right_slice;
		}

		// Equality relationships.
		template <index_t LhsHeight, index_t LhsWidth, typename LhsT,
				 index_t RhsHeight, index_t RhsWidth, typename RhsT>
		friend bool operator==(const matrix<LhsHeight, LhsWidth, LhsT>& lhs, 
				const matrix<RhsHeight, RhsWidth, RhsT>& rhs) noexcept
		{
			using commonT = std::common_type_t<LhsT, RhsT>;
			const commonT tolerance{equality_tolerance};
			
			if (LhsHeight != RhsHeight || LhsWidth != RhsWidth)
				return false;
			for (index_t r = 0; r < LhsHeight; ++r)
				for (index_t c = 0; c < LhsWidth; ++c)
				{
					if (std::abs(lhs[r][c] - rhs[r][c]) > tolerance)
						return false;
					if ((rhs[r][c] - lhs[r][c]) > tolerance)
						return false;
				}
			return true;
		}
		template <index_t LhsHeight, index_t LhsWidth, typename LhsT, 
				 index_t RhsHeight, index_t RhsWidth, typename RhsT>
		friend bool operator!=(const matrix<LhsHeight, LhsWidth, LhsT>& lhs, 
				const matrix<RhsHeight, RhsWidth, RhsT>& rhs) noexcept
		{
			return !(lhs == rhs);
		}

		// Matrix multiplication.
		template <index_t RhsWidth, typename RhsT>
		auto operator* (const matrix<Width, RhsWidth, RhsT>& rhs) noexcept
			-> matrix<Height, RhsWidth, std::common_type_t<T, RhsT>>
		{
			using commonT = std::common_type_t<T, RhsT>;
			constexpr index_t RhsHeight = Width; // Avoid confusion.
			matrix<RhsHeight, RhsWidth, commonT> product;
			
			for (index_t r = 0; r < RhsHeight; ++r)
				for (index_t c = 0; c < RhsWidth; ++c)
					product[r][c] = std::inner_product(
						storage[r].begin(), storage[r].end(), rhs.column_cbegin(c),	commonT(0)
					);
			return product;
		}
		
		// Streaming (printing).
		friend std::ostream& operator<<(std::ostream& stream, const self_t& matrix)
		{
			for (index_t r = 0; r < Height; ++r)
			{
				stream << '\n';
				for (index_t c = 0; c < Width; ++c)
					stream << '\t' << matrix[r][c];
			}
			return stream;
		}

		// Obtain an identity matrix. Only valid for square matrices.
		static self_t get_identity_matrix()
    	{
        	static_assert(Height == Width, "Identity matrix only defined for square matrices.");
			const T Size = Height;
			self_t identity {};
        	for (index_t i = 0; i < Size; ++i)
            	identity[i][i] = T(1);
        	return identity;
    	}
		
		// In place inversion. Only valid for square matrices.
		void invert()
		{
            static_assert(Height == Width, "Can only invert square matrices.");
			auto augmented_matrix = horizontal_concat(*this, get_identity_matrix());
            augmented_matrix.row_reduce();
            *this = augmented_matrix.get_right_slice();
        }

		// Inversion of the present matrix, returned by value.
		self_t get_inverse() const
		{
			self_t mtx{*this};
			mtx.invert();
			return mtx;
		}

		// Iteration through columns. Only the minimum feature set is implemented.
		class const_column_iterator 
		{
			const self_t& parent;
			index_t row;
			const index_t column;
			
			public:
			const_column_iterator(const self_t& parent, index_t row, index_t column) 
				: parent(parent), row(row), column(column) 
			{ }
			
			// Dereference.
			constexpr T operator*() const noexcept { return parent[row][column]; }
			
			// Increment.
			const_column_iterator& operator++() noexcept { ++row; return *this; }
			const_column_iterator operator++(int) noexcept
			{
				const_column_iterator unincremented{parent, row, column};
				++row;
				return unincremented;
			}

			// Decrement.
			const_column_iterator& operator--() noexcept { --row; return *this; }
			const_column_iterator operator--(int) noexcept
			{
				const_column_iterator undecremented{parent, row, column};
				--row;
				return undecremented;
			}
		}; // End of const_column_iterator.
	}; // End of class matrix.

	// Helper alias.
	template <index_t Size, typename T = default_T>
    using square_matrix = matrix<Size, Size, T>;

	// Concatenation, for producing the adjunct matrix.
    template <index_t Height, index_t LhsWidth, index_t RhsWidth, typename T = default_T>
    auto horizontal_concat(const matrix<Height, LhsWidth, T>& lhs, 
		const matrix<Height, RhsWidth, T>& rhs)
        -> matrix<Height, LhsWidth+RhsWidth, T>
    {
        matrix<Height, LhsWidth+RhsWidth, T> concatenation {};
        for (index_t r = 0; r < Height; ++r)
        {
            for (index_t c = 0; c < LhsWidth; ++c)
                concatenation[r][c] = lhs[r][c];
            for (index_t c = 0; c < RhsWidth; ++c)
                concatenation[r][c+LhsWidth] = rhs[r][c];
        }
        return concatenation; 
    }

} // End namespace matrix_math.

#endif // End ifndef CROWSTON_MATRIX_MATH_H.

