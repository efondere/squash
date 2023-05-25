/**
 * @file MatrixD.hpp
 * @author Eliot Fondere
 * @brief Dynamic Matrix class for mathematical calculations
 *
 * @copyright Copyright (c) 2023 Eliot Fondere (MIT License)
 */

#ifndef INCLUDE_SQH_MATRIX_D_HPP
#define INCLUDE_SQH_MATRIX_D_HPP

#include <cstddef>
#include <algorithm>
#include <functional>
#include <iostream>
#include <iomanip>

namespace sqh::math
{

template<typename T>
class MatrixD
{
public:
	T* data;

	size_t rowCount;
	size_t columnCount;

public:
	MatrixD(size_t rows, size_t columns)
	: rowCount(rows)
	, columnCount(columns)
	{
		data = new T[rows * columns];
	}

	MatrixD(const MatrixD<T>& other) = delete;
	MatrixD(MatrixD<T>& other) = delete;

	MatrixD(MatrixD<T>&& other) noexcept
    : data(other.data)
	, rowCount(other.rowCount)
	, columnCount(other.columnCount)
	{
		other.data = nullptr;
	}

	MatrixD(const MatrixD<T>&& other) = delete;
	MatrixD& operator=(const MatrixD<T>& other) = delete;
	MatrixD& operator=(MatrixD<T>& other) = delete;
	MatrixD& operator=(MatrixD<T>&& other) = delete;
	MatrixD& operator=(const MatrixD<T>&& other) = delete;

	MatrixD<double> operator-(const MatrixD<T>& other)
	{
		MatrixD<double> out(rowCount, columnCount);

		for (int i = 0; i < rowCount; i++)
		{
			for (int j = 0; j < columnCount; j++)
			{
				out.data[(columnCount * i) + j] = static_cast<double>(at(i, j)) - static_cast<double>(other.at(i, j));
			}
		}

		return out;
	}

	~MatrixD()
	{
		delete[] data;
	}

	void fromArray(T* array)
	{
		for (int i = 0; i < rowCount * columnCount; i++)
		{
			data[i] = array[i];
		}
	}

	inline T at(size_t i, size_t j) const
	{
		return data[(columnCount * i) + j];
	}

	double norm()
	{
		return sqrt(dot(*this));
	}

	double angle(MatrixD<T>& other)
	{
		return acos(dot(other) / (norm() * other.norm()));
	}

	double dot(MatrixD<T>& other)
	{
		double sum = 0.f;

		for (int i = 0; i < rowCount; i++)
		{
			for (int j = 0; j < columnCount; j++)
			{
				sum += static_cast<double>(at(i, j)) * static_cast<double>(other.at(i, j));
			}
		}

		return sum;
	}
};

} // sqh::math

#endif // INCLUDE_SQH_MATRIX_D_HPP
