/*
 * Matrix inversion benchmark (without atomics).
 *
 * Author: Robert H. Crowston, 2017.
 *
 *
 * Invoke with c++ -std=c++14 -O3
 *
 */

#include <chrono>
#include <iostream>
#include <random>

#include "matrix_math.hpp"

using counter_t = unsigned;
using timer = std::chrono::duration<double>;

using namespace matrix_math;

//
// time_random_matrices<>().
//
// Produces test_count random square matrices of dimension Size and tries to invert them.
// The total time spent in inversion functions is returned.
//
template <index_t Size>
timer time_random_matrices(unsigned test_count, counter_t& nonsingular_count, counter_t& degenerate_count)
{
	// Set up our source of random numbers. One seed per thread.
	std::random_device seed{};
	std::mt19937_64 generator{seed()};
	std::uniform_int_distribution<> distribution(-10, 10);

	timer inversion_time_elapsed{0};

	while (test_count-- > 0)
	{
		// Populate the next matrix with random elements.
		square_matrix<Size> matrix;
		for (auto& row : matrix)
			for (auto& element : row)
				element = distribution(generator);

		// Invert the matrix, if possible.
		auto start = std::chrono::high_resolution_clock::now();
		try
		{
			matrix.invert();
			++nonsingular_count;
		}
		catch (matrix_is_degenerate_error& )
		{
			++degenerate_count;
		}
		auto end = std::chrono::high_resolution_clock::now();
		inversion_time_elapsed += (end-start);
	}
	return inversion_time_elapsed;
}

int main ()
{
	// How many tests to run.
	const unsigned test_count = 1'000'000;

	// Counting.
	counter_t nonsingular_count {0};
	counter_t degenerate_count {0};
	
	timer inversion_time = time_random_matrices<7>(test_count, nonsingular_count, degenerate_count);

	std::cout << "Singular: " << nonsingular_count << "; " <<
		"degenerate: " << degenerate_count << ".\n" <<
		"Time spent in inversion functions: " << inversion_time.count() << " s.\n" <<
		"Average inversion time per matrix: " << (inversion_time.count() / test_count) << " s.\n";
	
	return 0;
}

