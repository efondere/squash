/**
 * @file squash_enc.hpp
 * @author Eliot Fondere
 * @brief
 *
 * @copyright Copyright (c) 2023 Eliot Fondere (MIT License)
 */

#ifndef INCLUDE_SQH_SQUASH_ENC_HPP
#define INCLUDE_SQH_SQUASH_ENC_HPP

#include "SquashHeader.hpp"
#include "../math/Matrix.hpp"
#include "../math/math.hpp"
#include "stb/stb_image.h"
#include <string>
#include <array>
#include <cmath>
#include <fstream>

namespace sqh
{

// TODO: this can go in a no-name namespace in a cpp file
namespace priv
{

// returns 1 if i == 0 and sqrt(2) otherwise
float delta_i(size_t i)
{
	if (i == 0) return 1.f;
	else return std::sqrt(2.f);
}

// gives the i, j entry of the C matrix for compression using DCT
float c_ij(size_t N, size_t i, size_t j)
{
	return (delta_i(i)/std::sqrt(static_cast<float>(N)))
	       * std::cos((static_cast<float>(i * (2 * j + 1)) * math::pi) / static_cast<float>(2 * N) );
}

// the horizontal and vertical dimensions of a block
constexpr size_t N = 8;

const math::Matrix<N, N, float> C =
	math::Matrix<N, N, float>::FromFunction([](size_t i, size_t j) {
		return c_ij(N, i, j);
	});

float c_haar_data[8][8] =
{
	{1.f/(2*sqrt(2.f)),  1.f/(2*sqrt(2.f)),  1.f/2,      0,  1/sqrt(2.f),             0,            0,            0},
	{1.f/(2*sqrt(2.f)),  1.f/(2*sqrt(2.f)),  1.f/2,      0, -1/sqrt(2.f),             0,            0,            0},
	{1.f/(2*sqrt(2.f)),  1.f/(2*sqrt(2.f)), -1.f/2,      0,             0,  1/sqrt(2.f),            0,            0},
	{1.f/(2*sqrt(2.f)),  1.f/(2*sqrt(2.f)), -1.f/2,      0,             0, -1/sqrt(2.f),            0,            0},
	{1.f/(2*sqrt(2.f)), -1.f/(2*sqrt(2.f)),      0,  1.f/2,             0,            0,  1/sqrt(2.f),            0},
	{1.f/(2*sqrt(2.f)), -1.f/(2*sqrt(2.f)),      0,  1.f/2,             0,            0, -1/sqrt(2.f),            0},
	{1.f/(2*sqrt(2.f)), -1.f/(2*sqrt(2.f)),      0, -1.f/2,             0,            0,            0,  1/sqrt(2.f)},
	{1.f/(2*sqrt(2.f)), -1.f/(2*sqrt(2.f)),      0, -1.f/2,             0,            0,            0, -1/sqrt(2.f)}
};

math::Matrix<N, N, float> c_haar = math::Matrix<N, N, float>::FromArray(c_haar_data);

// quantization table
uint8_t q_data[8][8] =
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
math::Matrix<N, N, float> q = math::Matrix<N, N, uint8_t>::FromArray(q_data).asType<float>();
math::Matrix<N, N, float> averageCoefficientsDCT{};
math::Matrix<N, N, float> averageCoefficientsHaar{};

} // namespace priv

// TODO: move this to a cpp file

void encode_haar(math::Matrix<8, 8, float> shifted_data)
{
	auto transformed_data = priv::c_haar.product(shifted_data.product(priv::c_haar.transpose()));
	priv::averageCoefficientsHaar = (priv::averageCoefficientsHaar + transformed_data.asType<float>([](float input) {return abs(input);})) * 0.5f;
}

std::array<int8_t, priv::N * priv::N> encode_block(math::Matrix<8, 8, uint8_t> block)
{

	// AKA f_tilde
	auto shifted_data = block.asType<float>([](uint8_t input) -> float {
		return static_cast<float>(input) - 128.f;
	});

	encode_haar(shifted_data);

	// AKA alpha_tilde (equation 12.9)
	auto transformed_data = priv::C.product(shifted_data.product(priv::C.transpose()));
	priv::averageCoefficientsDCT = (priv::averageCoefficientsDCT + transformed_data.asType<float>([](float input) {return abs(input);})) * 0.5f;

	auto quantized_data = (transformed_data / priv::q) + 0.5f;
	auto quantized_data_int = quantized_data.asType<int8_t>([](float input) {
		return static_cast<int8_t>(floorf(input));
	});

	return quantized_data_int.flatten([](size_t index)->math::MatrixIndex{
		return {index / priv::N, index % priv::N};
	});
}

void compress_block(const std::array<int8_t, priv::N * priv::N>& block_data, CompressedBlock& compressed_block)
{
	compressed_block.table = 0;
	compressed_block.dataCount = 0;

	for (int8_t value : block_data)
	{
		// shifting the first doesn't do anything
		compressed_block.table <<= 1;

		if (value != 0)
		{
			compressed_block.data[compressed_block.dataCount++] = value;
			compressed_block.table |= 1; // append a one at the end of the table, representing the value
		}
	}
}

int encode(const std::string& input_file_path, const std::string& output_file_path)
{
	int x, y, n;
	uint8_t* image_data = stbi_load(input_file_path.c_str(), &x, &y, &n, 3);

	auto output_file = std::ofstream(output_file_path, std::ios::binary);

	auto x_div = std::div(x, 8);
	auto y_div = std::div(y, 8);
	uint32_t x_blocks = x_div.quot + (x_div.rem == 0 ? 0 : 1);
	uint32_t y_blocks = y_div.quot + (y_div.rem == 0 ? 0 : 1);

	SquashHeader header{static_cast<uint32_t>(x), static_cast<uint32_t>(y), ImageChannels::RGB};
	output_file.write(reinterpret_cast<const char*>(&header), sizeof(header));
	output_file.write(reinterpret_cast<const char*>(priv::q_data), sizeof(uint8_t) * 8 * 8);
	//auto blockCount = x_blocks * y_blocks;
	//output_file.write(reinterpret_cast<const char*>(&blockCount), sizeof(blockCount));

	CompressedBlock compressed_block{0, {}, 0};

	for (int i = 0; i < x_blocks; i++) {
		for (int j = 0; j < y_blocks; j++) {
			for (int c = 0; c < 3; c++) {
				auto block = math::Matrix<priv::N, priv::N, uint8_t>::FromFunction(
						[&image_data, i, j, x, y, c](size_t k, size_t l) -> uint8_t {
							if (i + k >= x) return 128;
							if (j + l >= y) return 128;

							// TODO: something's suspicious here...
							return image_data[3 * ((x * ((priv::N * j) + l)) + (priv::N * i) + k) + c];
						});

				auto block_array = encode_block(block);
				compress_block(block_array, compressed_block);

				output_file.write(reinterpret_cast<const char *>(
										  &compressed_block.table), sizeof(compressed_block.table));
				output_file.write(reinterpret_cast<const char *>(
										  compressed_block.data), sizeof(int8_t) * compressed_block.dataCount);
			}
		}
	}

	output_file.close();
	stbi_image_free(image_data);

	priv::averageCoefficientsDCT.print();
	priv::averageCoefficientsHaar.print();

	return 0;
}

} // namespace sqh

#endif // INCLUDE_SQH_SQUASH_ENC_HPP
