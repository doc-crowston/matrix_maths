#include <iostream>

#include "matrix_math.hpp"

using namespace matrix_math;

int main ()
{
	try
	{
		square_matrix<3> mtx{ {1, 2, 3}, {0, 1, 0}, {5, 6, 0} };

		std::cout << "Matrix is " << mtx << std::endl;
		
		square_matrix<3> mtx2{mtx};
		
		mtx.invert();
		std::cout << "Inversion is " << mtx << std::endl;

		std::cout << "Multiplied is " << mtx2 * mtx << "\nor " << mtx * mtx2 << std::endl;

		return 0;
	}
	catch (std::exception& error)
	{
		std::cout << error.what() << std::endl;
	}

	/*
	 * To do.
	 *
	 *  2.  Add Catch()
	 *  3.  Unit test multiple 1x1, 2x2, 3x3, 4x4, 5x5, 6x6, 7x7 matrices (doubles).
	 *  4.  Improve code structure / fine tune implementation.
	 *  5.  Random number generator.
	 *  6.  Threading.
	 *  7.  Test performance on a few systems.
	 *  8.  Graphing.
	 *
	 */
}

