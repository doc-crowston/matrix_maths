/*
 * Matrix inversion benchmark.
 *
 * Author: Robert H. Crowston, 2017.
 *
 *
 * Invoke with c++ -std=c++14 -O3
 *
 */

#include <atomic>
#include <chrono>
#include <iostream>
#include <mutex>
#include <random>
#include <thread>
#include <vector>

#include "matrix_math.hpp"

using counter_t = std::atomic<unsigned>;
using timer = std::chrono::duration<double>;

using namespace matrix_math;

//
// time_random_matrices<>().
//
// Produces test_count random square matrices of dimension Size and tries to invert them.
// The total time spent in inversion functions is returned.
//
template <index_t Size>
timer time_random_matrices(unsigned test_count, counter_t& singular_count, counter_t& degenerate_count)
{
	// Set up our source of random numbers. One seed per thread.
	thread_local std::random_device seed{};
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
			++singular_count;
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
	const unsigned test_count = 50'000'000;
	// How many threads on which to execute these tests.	
	const unsigned thread_count = std::thread::hardware_concurrency();

	// Counting.
	counter_t singular_count {0};
	counter_t degenerate_count {0};
	timer inversion_time {0};
	std::mutex timer_mutex;		// Have to use a mutex since atomic<> won't support durations.

	// Thread management.
	std::vector<std::thread> pool;
	// Start threads.
	for (unsigned t = 0; t < thread_count; ++t)
	{
		pool.emplace_back( std::thread{ [&]
		{ 
			auto time = time_random_matrices<7>(
				test_count/thread_count, singular_count, degenerate_count
			);
			std::lock_guard<std::mutex> lock(timer_mutex);
			inversion_time += time;
		} } );
	}

	// Wait for threads to finish.
	for (auto& t : pool)
		t.join();

	std::cout << "Singular: " << singular_count << "; " <<
		"degenerate: " << degenerate_count << ".\n" <<
		"Time spent in inversion functions: " << inversion_time.count() << " s.\n" <<
		"Average inversion time per matrix: " << (inversion_time.count() / test_count) << " s.\n";
	
	return 0;
}

