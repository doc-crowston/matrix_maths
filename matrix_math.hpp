#include <array>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <ostream>
#include <numeric>
#include <type_traits>
#include <utility>

namespace matrix_math
{
    using index_t = std::size_t;
    using default_T = double;

    template <index_t Width, typename T = default_T>
    class row
    {
        using storage_t = std::array<T, Width>;
        storage_t storage{};

        public:
        using type = T;
        using self_type = row<Width, T>;
        
        row() noexcept : storage{} { }
        row(const std::initializer_list<T> init) noexcept 
		{
			index_t c = 0;
			for (const auto& init_col : init)
				storage[c++] = T(init_col);
		}
        
        constexpr T& operator[] (const index_t x) noexcept { return storage[x]; }
        constexpr const T operator[] (const index_t x) const noexcept { return storage[x]; }

		constexpr auto begin() noexcept { return storage.begin(); }
		constexpr auto end() noexcept { return storage.end(); }

        template <typename Callable>
        void apply(Callable f) noexcept(noexcept(f))
        // Apply a function f(element) to every element.
        {
            for (auto& element : storage)
                f(element);
        }

		//
		// Elementary row operations.
		// 

        // Row multiplication by a constant.
        void operator*= (const T rhs) { apply([rhs] (auto& el) {el *= rhs;}); }
        self_type operator* (const T rhs) const noexcept
        {
            self_type new_row{};
            for (index_t i = 0; i < Width; ++i)
                new_row[i] = storage[i] * rhs;
            return new_row;
        }

        // Addition of two rows.
        void operator+= (const row<Width,T>& rhs) noexcept
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
		using self_t = matrix<Height, Width, T>;

        private:
        using storage_t = std::array<row<Width, T>, Height>;
        storage_t storage{};

        public:
        matrix() noexcept : storage{} { }
        matrix(const std::initializer_list<std::initializer_list<T>> init) noexcept
		{
			index_t r = 0;
			for (const auto& row_init : init)
				storage[r++] = row_t{row_init};
		}

        constexpr row_t& operator[] (const index_t y) noexcept { return storage[y]; }
        constexpr const row_t& operator[] (const index_t y) const noexcept 
		{ 
			return storage[y]; 
		}

		constexpr auto begin() noexcept { return storage.begin(); }
		constexpr auto end() noexcept   { return storage.end(); }
		constexpr auto column_begin(index_t col) const noexcept 
		{ 
			return const_column_iterator(*this, 0, col);
		}
		constexpr auto column_end(index_t col) const noexcept
		{
			return const_column_iterator(*this, Height, col);
		}

		//
		// Elementary row operations.
		//
        void swap_rows(const index_t a, const index_t b) noexcept
        {
            using std::swap;
            swap(storage[a], storage[b]);
        }

		//
		// More complex operations.
		//
		void row_reduce()
		{
			// Following the Gauss--Jordan method.
			for (index_t r = 0; r < Height; ++r)
			{
				// We need a non-zero element at position [r][r].
				// Interchange rows to obtain a non-zero.
				if (storage[r][r] == T(0))
				{
					index_t s = r+1;
					while (s < Height && storage[s][s] == T(0))
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
		template <index_t RhsHeight, index_t RhsWidth, typename RhsT>
		bool operator==(const matrix<RhsHeight, RhsWidth, RhsT>& rhs) noexcept
		{
			if (RhsHeight != Height || RhsWidth != Width)
				return false;
			for (index_t r = 0; r < Height; ++r)
				for (index_t c = 0; c < Width; ++c)
					if (storage[r][c] != rhs[r][c])
							return false;
			return true;
		}
		template <index_t RhsHeight, index_t RhsWidth, typename RhsT>
		bool operator!=(const matrix<RhsHeight, RhsWidth, RhsT>& rhs) noexcept
		{
			return !(*this == rhs);
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
					product[r][c] = std::inner_product(storage[r].begin(), storage[r].end(), rhs.column_begin(c), commonT(0));
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

		//
		// Square matrix operations.
		//
		static self_t get_identity_matrix()
    	{
        	static_assert(Height == Width);
			const T Size = Height;
			self_t identity {};
        	for (index_t i = 0; i < Size; ++i)
            	identity[i][i] = T(1);
        	return identity;
    	}

		void invert()
		{
            static_assert(Height == Width);
			const T Size = Height; // == Width
			auto augmented_matrix = horizontal_concat(*this, get_identity_matrix());
            augmented_matrix.row_reduce();
            *this = augmented_matrix.get_right_slice();
        }

		//
		// Iteration through columns.
		//
		class const_column_iterator 
		{
			const self_t& parent;
			index_t row;
			const index_t column;
			
			public:
			const_column_iterator(const self_t& parent, index_t row, index_t column) : parent(parent), row(row), column(column) { }
			
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
		};

		//
		// Inversion fail.
		//
		class matrix_is_degenerate_error : public std::domain_error
		{
			public:
			matrix_is_degenerate_error() : std::domain_error("Cannot invert degenerate matrix.") {}
			virtual ~matrix_is_degenerate_error() {}
		};
    }; 

	template <index_t Size, typename T = default_T>
    using square_matrix = matrix<Size, Size, T>;
 
    template <index_t Height, index_t LhsWidth, index_t RhsWidth, typename T = default_T>
    auto horizontal_concat(const matrix<Height, LhsWidth, T>& lhs, const matrix<Height, RhsWidth, T>& rhs)
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

}

