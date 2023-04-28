/**
 * @file main.cpp
 * @author Eliot Fondere
 * @brief
 *
 * @copyright Copyright (c) 2023 Eliot Fondere (MIT License)
 */

#include "math/Matrix.hpp"
#include "math/math.hpp"
#include <iostream>
#include <cmath>

namespace mat = sqh::math;

float delta_i(size_t i)
{
	if (i == 0) return 1.f;
	else return std::sqrt(2.f);
}

float c_ij(size_t N, size_t i, size_t j)
{
	return (delta_i(i)/std::sqrt(static_cast<float>(N)))
		* std::cos((static_cast<float>(i * (2 * j + 1)) * mat::pi) / static_cast<float>(2 * N) );
}

int main()
{
	const size_t N = 8;

    mat::Matrix<N, N, float> m =
			mat::Matrix<N, N, float>::FromFunction([](size_t i, size_t j){return c_ij(N, i, j);});
	m.print();

	return 0;
}
