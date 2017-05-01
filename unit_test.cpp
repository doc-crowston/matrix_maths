#include "matrix_math.hpp"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace matrix_math;

TEST_CASE( "1x1 matrices.", "[1x1]" ) 
{

	SECTION( "The identity matrix." )
	{
		square_matrix<1> mtx{ {1} };
		REQUIRE( mtx == mtx.get_inverse() );
	}
	SECTION( "Arbitrary matrix 1." )
	{
		square_matrix<1> mtx{ {103217.4} };

		REQUIRE( mtx * mtx.get_inverse() == mtx.get_inverse() * mtx );
		REQUIRE( mtx.get_inverse() * mtx == mtx * mtx.get_inverse() );
		REQUIRE( mtx.get_inverse() * mtx == square_matrix<1>::get_identity_matrix() );
	}
	SECTION( "Arbitrary matrix 2." )
	{
		square_matrix<1> mtx{ {-103213217.4} };

		REQUIRE( mtx * mtx.get_inverse() == mtx.get_inverse() * mtx );
		REQUIRE( mtx.get_inverse() * mtx == square_matrix<1>::get_identity_matrix() );
	}
	SECTION( "Arbitrary matrix 3." )
	{
		square_matrix<1> mtx{ {0.004124} };

		REQUIRE( mtx * mtx.get_inverse() == mtx.get_inverse() * mtx );
		REQUIRE( mtx.get_inverse() * mtx == square_matrix<1>::get_identity_matrix() );
	}
}

TEST_CASE( "2x2 matrices.", "[2x2]" )
{
	const auto identity = square_matrix<2>::get_identity_matrix();

	SECTION( "Invertible matrices." )
	{
		REQUIRE ( identity == identity.get_inverse() );
		
		square_matrix<2> mtx_1{ {2, 7}, {4, 6} };
		REQUIRE( mtx_1.get_inverse() * mtx_1 == identity );
	
		square_matrix<2> mtx_2{ {0, 1}, {1, 2} };
		REQUIRE( mtx_2.get_inverse() * mtx_2 == identity );

		square_matrix<2> mtx_3{ {0.7, 1.99}, {24.1, 9999} };
		REQUIRE( mtx_3.get_inverse() * mtx_3 == identity );
	}
	
	SECTION( "Degenerate matrices." )
	{
		square_matrix<2> mtx_4{ {0, 1}, {1, 0} };
		CHECK_THROWS(mtx_4.invert());
		
		square_matrix<2> mtx_5{ {2, 6}, {1, 3} };
		CHECK_THROWS(mtx_5.invert());
		
		square_matrix<2> mtx_6{ {10, 10}, {10, 10} };
		CHECK_THROWS(mtx_6.invert());
		
		square_matrix<2> mtx_7{ {0.001, 0.002}, {0.003, 0.006} };
		CHECK_THROWS(mtx_7.invert());
	}
}

TEST_CASE( "3x3 matrices.", "[3x3]" )
{
	const auto identity = square_matrix<3>::get_identity_matrix();

	SECTION( "Invertible matrices." )
	{
		{
			square_matrix<3> mtx{ 
				{-1,  3, -3}, 
				{ 0, -6,  5},
				{-5, -3,  1}
			};
			REQUIRE( mtx.get_inverse() * mtx == identity );
		}
		{
			square_matrix<3> mtx{ 
				{ 7,  2,  1}, 
				{ 0,  3, -1},
				{-3,  4, -2}
			};
			REQUIRE( mtx.get_inverse() * mtx == identity );
		}
		{
			square_matrix<3> mtx{ 
				{ 2,  1,  0}, 
				{ 0,  2,  0},
				{ 2,  0,  1}
			};
			REQUIRE( mtx.get_inverse() * mtx == identity );
		}
	}

	SECTION( "Degenerate matrices." )
	{
		{
			square_matrix<3> mtx{ 
				{ 1,  0,  0}, 
				{-2,  0,  0},
				{ 4,  6,  1}
			};
			CHECK_THROWS( mtx.get_inverse() * mtx == identity );
		}
		{
			square_matrix<3> mtx{ 
				{ 0,  1,  0}, 
				{ 0,  0,  1},
				{ 0,  1,  0}
			};
			CHECK_THROWS( mtx.get_inverse() * mtx == identity );
		}
	}	

}

TEST_CASE( "4x4 matrices.", "[4x4]" )
{
	SECTION( "Invertible matrices." )
	{
		const auto identity = square_matrix<4>::get_identity_matrix();
		{
			square_matrix<4> mtx{
				{ 4,  0,  0,  0},
				{ 0,  0,  2,  0},
				{ 0,  1,  2,  0},
				{ 1,  0,  0,  1}
			};
			REQUIRE( mtx * mtx.get_inverse() == identity );
		}
		{
			square_matrix<4> mtx{
				{ 1,  2,  1,  0},
				{ 2,  1,  1,  1},
				{-1,  2,  1, -1},
				{ 1,  1,  1,  2}
			};
			REQUIRE( mtx * mtx.get_inverse() == identity );
		}	
	}
}

TEST_CASE( "7x7 matrices.", "[7x7]" )
{
	SECTION( "Invertible matrices." )
	{
		const auto identity = square_matrix<7>::get_identity_matrix();
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
			REQUIRE( mtx * mtx.get_inverse() == identity );
		}
	}
}

