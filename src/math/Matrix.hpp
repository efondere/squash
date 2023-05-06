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

template<typename T, typename N_T>
inline N_T default_conversion(T input)
{
	return static_cast<N_T>(input);
}

template<size_t R, size_t C, typename T>
class Matrix
{
public:
	struct Index
	{
		size_t i;
		size_t j;
	};

	Matrix();
	static Matrix<R, C, T> FromFunction(T(*function)(size_t, size_t));
	static Matrix<R, C, T> FromArray(T array[R][C]);
	static Matrix<R, C, T> Outer(Matrix<1, C, T> rows, Matrix<R, 1, T> columns);
	static Matrix<R, C, T> Identity();
	// static Transpose() -> or have it as a separate function?

	Matrix<R, 1, T> getColumn(size_t index);
	Matrix<1, C, T> getRow(size_t index);
	Matrix<C, 1, T> getRowT(size_t index);

	std::array<T, R * C>&& flatten(Index(*flatten_function)(size_t));

	template<typename N_T>
	Matrix<C, R, N_T> asType(N_T(*conversion_function)(T input)=default_conversion);

	Matrix<C, R, T> transpose();

	Matrix<R, C, T> operator-(T scalar);
	Matrix<R, C, T> operator*(T scalar);
	Matrix<R, C, T> operator/(Matrix<R, C, T>& other);

	template<size_t C_2>
	Matrix<R, C_2, T> product(Matrix<C, C_2, T> other);

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

template<size_t R, size_t C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::Outer(Matrix<1, C, T> rows, Matrix<R, 1, T> columns)
{
	Matrix<R, C, T> m;

	for (size_t i = 0; i < R; i++)
	{
		for (size_t j = 0; j < C; j++)
		{
			m.data[i][j] = columns.data[i][0] * rows.data[0][j];
		}
	}

	return m;
}

template<size_t R, size_t C, typename T>
Matrix<R, 1, T> Matrix<R, C, T>::getColumn(size_t index)
{
	Matrix<R, 1, T> m;

	for (size_t i = 0; i < R; i++)
	{
		m.data[i][0] = data[i][index];
	}

	return m;
}

template<size_t R, size_t C, typename T>
Matrix<1, C, T> Matrix<R, C, T>::getRow(size_t index)
{
	Matrix<1, C, T> m;

	for (size_t i = 0; i < R; i++)
	{
		m.data[0][i] = data[index][i];
	}

	return m;
}

template<size_t R, size_t C, typename T>
Matrix<C, 1, T> Matrix<R, C, T>::getRowT(size_t index)
{
	Matrix<1, C, T> m;

	for (size_t i = 0; i < R; i++)
	{
		m.data[i][0] = data[index][i];
	}

	return m;
}

template<size_t R, size_t C, typename T>
std::array<T, R * C>&& Matrix<R, C, T>::flatten(Matrix::Index (*flatten_function)(size_t)) {
	static_assert(false, "flatten no implemented");
	return std::array<T, R * C>();
}

template<size_t R, size_t C, typename T>
Matrix<C, R, T> Matrix<R, C, T>::transpose() {
	Matrix<C, R, T> m;

	for (size_t i = 0; i < C; i++)
	{
		for (size_t j = 0; j < R; j++)
		{
			m.data[i][j] = data[j][i];
		}
	}

	return m;
}

template<size_t R, size_t C, typename T>
template<size_t C_2>
Matrix<R, C_2, T> Matrix<R, C, T>::product(Matrix<C, C_2, T> other)
{
	Matrix<R, C_2, T> m;

	for (size_t i = 0; i < R; i++)
	{
		for (size_t j = 0; j < C_2; j++)
		{
			T sum = 0;
			for (size_t k = 0; k < C; k++)
			{
				sum += data[i][k] * other.data[k][j];
			}

			m.data[i][j] = sum;
		}
	}

	return m;
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::FromArray(T array[R][C]) {
	Matrix<R, C, T> m;

	for (size_t j = 0; j < C; j++)
	{
		for (size_t i = 0; i < R; i++)
		{
			m.data[i][j] = array[i][j];
		}
	}

	return m;
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::operator*(T scalar)
{
	Matrix<R, C, T> m;

	for (size_t j = 0; j < C; j++)
	{
		for (size_t i = 0; i < R; i++)
		{
			m.data[i][j] = data[i][j] * scalar;
		}
	}

	return m;
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::operator-(T scalar)
{
	Matrix<R, C, T> m;

	for (size_t j = 0; j < C; j++)
	{
		for (size_t i = 0; i < R; i++)
		{
			m.data[i][j] = data[i][j] - scalar;
		}
	}

	return m;
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::operator/(Matrix<R, C, T>& other) {
	Matrix<R, C, T> m;

	for (size_t j = 0; j < C; j++)
	{
		for (size_t i = 0; i < R; i++)
		{
			m.data[i][j] = data[i][j] / other.data[i][j];
		}
	}

	return m;
}

template<size_t R, size_t C, typename T>
template<typename N_T>
Matrix<C, R, N_T> Matrix<R, C, T>::asType(N_T (*conversion_function)(T))
{
	Matrix<C, R, N_T> m;

	for (size_t j = 0; j < C; j++)
	{
		for (size_t i = 0; i < R; i++)
		{
			m.data[i][j] = conversion_function(data[i][j]);
		}
	}

	return m;
}

template<size_t R, size_t C, typename T>
constexpr Matrix<R, C, T> operator*(T scalar, Matrix<R, C, T> matrix)
{
	Matrix<R, C, T> m;

	for (size_t j = 0; j < C; j++)
	{
		for (size_t i = 0; i < R; i++)
		{
			m.data[i][j] = matrix.data[i][j] * scalar;
		}
	}

	return m;
}

} // math
} // sqh

#endif // INCLUDE_SQH_MATRIX_HPP
