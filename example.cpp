/*
 * Matrix inversion: example of API.
 *
 * Author: Robert H. Crowston, 2017.
 *
 *
 * Invoke with c++ -std=c++14 -O3
 *
 */

#include <iomanip>
#include <iostream>

#include "matrix_math.hpp"

int main ()
{
	using namespace matrix_math;
	try
	{
		// Create and populate a 7x7 square matrix.
		square_matrix<7> mtx{
			{ 1,  2,  3,  4,  0, -1,  0},
			{ 0,  1,  1,  0,  1,  0,  0},
			{ 1,  0,  0,  0,  0,  1,  0},
			{ 0,  2,  2,  2, -2,  1,  3},
			{ 1,  3,  5,  7,  0, -1,  1},
			{ 0,  0,  1,  0,  1,  0,  0},
			{ 9, -2,  0,  0,  0,  2,  0}
		};
		
		// Native support for streaming the matrix is provided.
		std::cout << "Matrix is " << mtx << std::endl;
		
		// Inversion is availble through the ::get_inverse() or the (in place) ::invert() methods.
		std::cout << "Inversion is " << std::fixed << mtx.get_inverse() << std::endl;

		std::cout << "Multiplied is " << mtx.get_inverse() * mtx << "\nor " << mtx * mtx.get_inverse() << std::endl;

		return 0;
	}
	catch (matrix_is_degenerate_error& error)
	{
		// If the matrix cannot be inverted, an exeception will be thrown.
		std::cout << "Cannot invert: " << error.what() << std::endl;
	}
	catch (std::exception& error)
	{
		std::cout << error.what() << std::endl;
	}
}

