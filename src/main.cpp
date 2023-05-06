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

int main() {
	const size_t N = 8;

	mat::Matrix<N, N, float> C =
			mat::Matrix<N, N, float>::FromFunction([](size_t i, size_t j) { return c_ij(N, i, j); });

	float f_data[8][8] =
		{
			{40, 193, 89, 37, 209, 236, 41, 14},
			{102, 165, 36, 150, 247, 104, 7, 19},
			{157, 92, 88, 251, 156, 3, 20, 35},
			{153, 75, 220, 193, 29, 13, 34, 22},
			{116, 173, 240, 54, 11, 38, 20, 19},
			{162, 255, 109, 9, 26, 22, 20, 29},
			{237, 182, 5, 28, 20, 15, 28, 20},
			{222, 33, 8, 23, 24, 29, 23, 23}
		};

	float q_data[8][8] =
		{
			{10, 16, 22, 28, 34, 40, 46, 52},
			{16, 22, 28, 34, 40, 46, 52, 58},
			{22, 28, 34, 40, 46, 52, 58, 64},
			{28, 34, 40, 46, 52, 58, 64, 70},
			{34, 40, 46, 52, 58, 64, 70, 76},
			{40, 46, 52, 58, 64, 70, 76, 82},
			{46, 52, 58, 64, 70, 76, 82, 88},
			{52, 58, 64, 70, 76, 82, 88, 94}
		};

	mat::Matrix<N, N, float> f = mat::Matrix<N, N, float>::FromArray(f_data);
	auto f_tilde = f - 128;
	auto alpha_tilde = C.product((f_tilde.product(C.transpose())));


	mat::Matrix<N, N, float> q = mat::Matrix<N, N, float>::FromArray(q_data);
	auto l = alpha_tilde / q - (-0.5);
	mat::Matrix<N, N, int> l_int = l.asType<int>([](float input){return static_cast<int>(floorf(input));});

	l_int.print();

	return 0;
}
