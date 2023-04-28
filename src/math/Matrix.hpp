/**
 * @file Matrix.hpp
 * @author Eliot Fondere
 * @brief Matrix class for mathematical calculations
 *
 * @copyright Copyright (c) 2023 Eliot Fondere (MIT License)
 */

#ifndef INCLUDE_SQH_MATRIX_HPP
#define INCLUDE_SQH_MATRIX_HPP

#include <cstddef>
#include <algorithm>

#include <iostream>
#include <iomanip>

namespace sqh
{
namespace math
{

template<size_t R, size_t C, typename T>
class Matrix
{
public:
	Matrix();
	static Matrix<R, C, T> FromFunction(T(*function)(size_t, size_t));
	static Matrix<R, C, T> Outer(Matrix<1, C, T> rows, Matrix<R, 1, T> columns);
	static Matrix<R, C, T> Identity();
	// static Transpose() -> or have it as a separate function?

	Matrix<R, 1, T> getColumn(size_t index);
	Matrix<1, C, T> getRow(size_t index);
	Matrix<C, 1, T> getRowT(size_t index);

	void print();

public:
	T data[R][C];
};

template<size_t R, size_t C, typename T>
void Matrix<R, C, T>::print()
{
	const int initial_precision = static_cast<int>(std::cout.precision());
	std::cout << std::right << std::fixed << std::setprecision(3);

	for (size_t i = 0; i < R; i++)
	{
		if (R == 1)
			std::cout << "[";
		else if (i == 0)
			std::cout << "/";
		else if (i == R - 1)
			std::cout << "\\";
		else
			std::cout << "|";

		for (size_t j = 0; j < C; j++)
		{
			std::cout << " " << std::setw(7) << data[i][j];
		}

		if (R == 1)
			std::cout << " ]";
		else if (i == 0)
			std::cout << " \\";
		else if (i == R - 1)
			std::cout << " /";
		else
			std::cout << " |";

		std::cout << std::endl;
	}

	std::cout << std::left << std::defaultfloat << std::setprecision(initial_precision);
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T>::Matrix()
: data{}
{
	static_assert(R != 0, "Matrix cannot have 0 rows.");
	static_assert(C != 0, "Matrix cannot have 0 columns.");
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::Identity()
{
	static_assert(R == C, "Matrix is not square. (At Matrix::Identity)");

	Matrix<R, C, T> m;
	for (size_t i = 0; i < R; i++)
	{
		m.data[i][i] = static_cast<T>(1);
	}

	return m;
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::FromFunction(T (*function)(size_t, size_t))
{
	Matrix<R, C, T> m;

	for (size_t i = 0; i < R; i++)
	{
		for (size_t j = 0; j < C; j++)
		{
			m.data[i][j] = function(i, j);
		}
	}

	return m;
}

} // math
} // sqh

#endif // INCLUDE_SQH_MATRIX_HPP
