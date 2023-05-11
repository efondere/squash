/**
 * @file SquashHeader.hpp
 * @author Eliot Fondere
 * @brief Header definition for all squash files
 *
 * @copyright Copyright (c) 2023 Eliot Fondere (MIT License)
 */

#ifndef INCLUDE_SQH_SQUASH_HEADER_HPP
#define INCLUDE_SQH_SQUASH_HEADER_HPP

#include <cstdint>

namespace sqh
{

enum class ImageChannels: uint8_t
{
	Grey,
	RGB,
	RGBA
};

struct SquashHeader
{
	uint32_t sizeX;
	uint32_t sizeY;

	ImageChannels channels;
};

// TODO: store magic number
struct CompressedBlock
{
	uint64_t table;
	int8_t data[64]; // block size (N) should be declared in this file and this should be N * N instead of 64
	uint8_t dataCount;
};

} // namespace sqh

#endif // INCLUDE_SQH_SQUASH_HEADER_HPP
