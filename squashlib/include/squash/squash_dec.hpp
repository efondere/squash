/**
 * @file squash_dec.hpp
 * @author Eliot Fondere
 * @brief
 *
 * @copyright Copyright (c) 2023 Eliot Fondere (MIT License)
 */

#ifndef INCLUDE_SQH_SQUASH_DEC_HPP
#define INCLUDE_SQH_SQUASH_DEC_HPP

#include "SquashHeader.hpp"
#include "stb/stb_image_write.h"
#include "../math/Matrix.hpp"
#include <string>
#include <fstream>
#include <iostream>

namespace sqh
{

namespace priv
{
    constexpr uint64_t mask = uint64_t(1) << 63;
}

void decompress_block(uint64_t table, std::ifstream& input_file, int8_t block_data[8][8])
{
    for (int k = 0; k < 64; k++)
    {
        auto index = std::div(k, 8);

        if (!(table & priv::mask))
        {
            // TODO: figure out where the matrix is being transposed... Check Matrix::FromArray
            // TODO: if mask[63] = 0 we can tag the block as being special
            block_data[index.rem][index.quot] = 0;
        }
        else
        {
            int8_t value;
            input_file.read(reinterpret_cast<char*>(&value), sizeof(value));
            block_data[index.rem][index.quot] = value;
        }

            table <<= 1;
    }
}

math::Matrix<8, 8, uint8_t> decode_block(std::ifstream& input_file, math::Matrix<8, 8, float>& q)
{
	uint64_t table;
	input_file.read(reinterpret_cast<char*>(&table), sizeof(table));

	int8_t block_data[8][8];
	decompress_block(table, input_file, block_data);

	auto block = math::Matrix<8, 8, int8_t>::FromArray(block_data);

	auto beta = block.asType<float>() * q;
	auto fbarf = priv::C.transpose().product(beta.product(priv::C)) + 128.f;
	auto fbar = fbarf.asType<uint8_t>([](float input) {
		return static_cast<uint8_t>(std::min(255.f, std::max(0.f, floorf(input))));
	});

	return fbar;
}

int decode(const std::string& input_file_path, const std::string& output_file_path)
{
	std::ifstream input_file(input_file_path, std::ios::binary);

	SquashHeader file_header{0, 0, ImageChannels::Grey};
	input_file.read(reinterpret_cast<char*>(&file_header), sizeof(file_header));
	std::cout << file_header.sizeX << " - " << file_header.sizeY << " - " << int(file_header.channels) << std::endl;
	uint8_t q_data[8][8] = {};
	input_file.read(reinterpret_cast<char*>(q_data), sizeof(uint8_t) * 8 * 8);
	auto q = math::Matrix<8, 8, uint8_t>::FromArray(q_data).asType<float>();
	q.asType<int>().print();

	auto x_div = std::div(file_header.sizeX, 8);
	auto y_div = std::div(file_header.sizeY, 8);
	uint32_t x_blocks = x_div.quot + (x_div.rem == 0 ? 0 : 1);
	uint32_t y_blocks = y_div.quot + (y_div.rem == 0 ? 0 : 1);

	std::cout << x_blocks << " - " << y_blocks << std::endl;

	auto* image_data = new uint8_t[file_header.sizeX * file_header.sizeY * 3]; // 3 is the amount of channels
	uint64_t mask = 1; mask <<= 63;

	for (int i = 0; i < x_blocks; i++) {
		for (int j = 0; j < y_blocks; j++) {
			for (int c = 0; c < 3; c++) {
				auto fbar = decode_block(input_file, q);

				for (int k = 0; k < 8; k++) {
					for (int l = 0; l < 8; l++) {
						if ((k + i * 8 >= file_header.sizeX) || (l + j * 8 >= file_header.sizeY)) continue;
						// TODO: why is it transposed?
						image_data[3* ((file_header.sizeX * ((priv::N * j) + l)) + (priv::N * i) + k) + c] = fbar.data[l][k];
					}
				}
			}
		}
	}

	stbi_write_png(output_file_path.c_str(), file_header.sizeX, file_header.sizeY, 3, &image_data[0], 0);
	delete[] image_data;
	input_file.close();
	return 0;
}

} // namespace sqh

#endif // INCLUDE_SQH_SQUASH_DEC_HPP
