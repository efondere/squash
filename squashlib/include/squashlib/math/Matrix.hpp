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
#include <functional>
#include <iostream>
#include <iomanip>

namespace sqh::math
{

template<typename T, typename N_T>
inline N_T default_conversion(T input)
{
	return static_cast<N_T>(input);
}

struct MatrixIndex
{
	size_t i;
	size_t j;
};

template<size_t R, size_t C, typename T>
class Matrix
{
public:
	Matrix();
	static Matrix<R, C, T> FromFunction(std::function<T(size_t, size_t)> function);
	static Matrix<R, C, T> FromArray(T array[R][C]);
	static Matrix<R, C, T> Outer(Matrix<1, C, T> rows, Matrix<R, 1, T> columns);
	static Matrix<R, C, T> Identity();
	static Matrix<R, C, T> Transpose(const Matrix<C, R, T>& other);

	Matrix<R, 1, T> getColumn(size_t index);
	Matrix<1, C, T> getRow(size_t index);
	Matrix<C, 1, T> getRowT(size_t index);

	std::array<T, R * C> flatten(const std::function<MatrixIndex(size_t)>& flatten_function) const;
	std::array<T, R * C> flatten(const std::function<MatrixIndex(size_t, size_t, size_t)>& flatten_function) const;

	template<typename N_T>
	Matrix<C, R, N_T> asType2(const std::function<N_T(T)>& conversion_function) const;

	template<typename N_T>
	Matrix<C, R, N_T> asType(N_T(*conversion_function)(T)=default_conversion) const;

	Matrix<C, R, T> transpose() const;

	Matrix<R, C, T> operator+(T scalar);
	Matrix<R, C, T> operator-(T scalar);
	Matrix<R, C, T> operator*(T scalar);
	Matrix<R, C, T> operator/(T scalar);

	Matrix<R, C, T> operator+(const Matrix<R, C, T>& other);
	Matrix<R, C, T> operator-(const Matrix<R, C, T>& other);
	Matrix<R, C, T> operator*(const Matrix<R, C, T>& other);
	Matrix<R, C, T> operator/(const Matrix<R, C, T>& other);

	template<size_t C_2>
	Matrix<R, C_2, T> product(const Matrix<C, C_2, T>& other) const;

	T dot(const Matrix<R, C, T>& other) const;
	T norm() const;
	T normalized_dot(const Matrix<R, C, T>& other) const;

	void print();
	void toStream(std::ostream& output_stream) const;

public:
	T data[R][C];
};


// IMPLEMENTATION
template<size_t R, size_t C, typename T>
Matrix<R, C, T>::Matrix()
	: data{}
{
	static_assert(R != 0, "Matrix cannot have 0 rows.");
	static_assert(C != 0, "Matrix cannot have 0 columns.");
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::FromFunction(std::function<T(size_t, size_t)> function)
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
Matrix<R, C, T> Matrix<R, C, T>::Transpose(const Matrix<C, R, T> &other)
{
	return other.transpose();
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
std::array<T, R * C> Matrix<R, C, T>::flatten(const std::function<MatrixIndex(size_t)>& flatten_function) const
{
	std::array<T, R * C> array;

	for (size_t i = 0; i < array.size(); i++)
	{
		MatrixIndex index = flatten_function(i);
		array.at(i) = data[index.i][index.j];
	}

	return array;
}

template<size_t R, size_t C, typename T>
std::array<T, R * C> Matrix<R, C, T>::flatten(const std::function<MatrixIndex(size_t, size_t, size_t)>& flatten_function) const
{
	std::array<T, R * C> array;

	for (size_t i = 0; i < array.size(); i++)
	{
		MatrixIndex index = flatten_function(i, R, C);
		array.at(i) = data[index.i][index.j];
	}

	return array;
}

template<size_t R, size_t C, typename T>
template<typename N_T>
Matrix<C, R, N_T> Matrix<R, C, T>::asType2(const std::function<N_T(T)>& conversion_function) const
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
template<typename N_T>
Matrix<C, R, N_T> Matrix<R, C, T>::asType(N_T(*conversion_function)(T)) const
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
Matrix<C, R, T> Matrix<R, C, T>::transpose() const {
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
Matrix<R, C, T> Matrix<R, C, T>::operator+(T scalar)
{
	Matrix<R, C, T> m;

	for (size_t j = 0; j < C; j++)
	{
		for (size_t i = 0; i < R; i++)
		{
			m.data[i][j] = data[i][j] + scalar;
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
Matrix<R, C, T> Matrix<R, C, T>::operator/(T scalar)
{
	Matrix<R, C, T> m;

	for (size_t j = 0; j < C; j++)
	{
		for (size_t i = 0; i < R; i++)
		{
			m.data[i][j] = data[i][j] / scalar;
		}
	}

	return m;
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::operator+(const Matrix<R, C, T>& other) {
	Matrix<R, C, T> m;

	for (size_t j = 0; j < C; j++)
	{
		for (size_t i = 0; i < R; i++)
		{
			m.data[i][j] = data[i][j] + other.data[i][j];
		}
	}

	return m;
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::operator-(const Matrix<R, C, T>& other) {
	Matrix<R, C, T> m;

	for (size_t j = 0; j < C; j++)
	{
		for (size_t i = 0; i < R; i++)
		{
			m.data[i][j] = data[i][j] - other.data[i][j];
		}
	}

	return m;
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::operator*(const Matrix<R, C, T>& other) {
	Matrix<R, C, T> m;

	for (size_t j = 0; j < C; j++)
	{
		for (size_t i = 0; i < R; i++)
		{
			m.data[i][j] = data[i][j] * other.data[i][j];
		}
	}

	return m;
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> Matrix<R, C, T>::operator/(const Matrix<R, C, T>& other) {
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
T Matrix<R, C, T>::dot(const Matrix<R, C, T>& other) const
{
	T sum = 0;
	for (int i = 0; i < R; i++)
	{
		for (int j = 0; j < C; j++)
		{
			sum += data[i][j] * other.data[i][j];
		}
	}

	return sum;
}

template<size_t R, size_t C, typename T>
T Matrix<R, C, T>::norm() const
{
	return sqrt(dot(*this));
}

template<size_t R, size_t C, typename T>
T Matrix<R, C, T>::normalized_dot(const Matrix<R, C, T> &other) const
{
	return dot(other) / (norm() * other.norm());
}

template<size_t R, size_t C, typename T>
template<size_t C_2>
Matrix<R, C_2, T> Matrix<R, C, T>::product(const Matrix<C, C_2, T>& other) const
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
			std::cout << " " << std::setw(8) << data[i][j];
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
void Matrix<R, C, T>::toStream(std::ostream &output_stream) const
{
	for (size_t i = 0; i < R; i++)
	{
		for (size_t j = 0; j < C; j++)
		{
			if (j != 0)
				output_stream << " ";

			output_stream << data[i][j];
		}

		output_stream << "\n";
	}

	output_stream << std::endl;
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

template<size_t R, size_t C, typename T>
Matrix<R, C, T> operator/(T scalar, const Matrix<R, C, T>& matrix)
{
	Matrix<R, C, T> m;

	for (size_t j = 0; j < C; j++)
	{
		for (size_t i = 0; i < R; i++)
		{
			m.data[i][j] = scalar / matrix.data[i][j];
		}
	}

	return m;
}

} // sqh::math

#endif // INCLUDE_SQH_MATRIX_HPP
