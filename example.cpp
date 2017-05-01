#include <iomanip>
#include <iostream>

#include "matrix_math.hpp"

using namespace matrix_math;

int main ()
{
	try
	{
		square_matrix<7> mtx{
			{ 1,  2,  3,  4,  0, -1,  0},
			{ 0,  1,  1,  0,  1,  0,  0},
			{ 1,  0,  0,  0,  0,  1,  0},
			{ 0,  2,  2,  2, -2,  1,  3},
			{ 1,  3,  5,  7,  0, -1,  1},
			{ 0,  0,  1,  0,  1,  0,  0},
			{ 9, -2,  0,  0,  0,  2,  0}
		};

		std::cout << "Matrix is " << mtx << std::endl;
		
		std::cout << "Inversion is " << std::fixed << mtx.get_inverse() << std::endl;

		std::cout << "Multiplied is " << mtx.get_inverse() * mtx << "\nor " << mtx * mtx.get_inverse() << std::endl;

		return 0;
	}
	catch (std::exception& error)
	{
		std::cout << error.what() << std::endl;
	}

}

