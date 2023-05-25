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

constexpr size_t BLOCK_SIZE = 8;
constexpr uint32_t MAGIC_NUMBER = 0x2F737168;

enum class ImageChannels: uint8_t
{
	Grey = 1,
	RGB  = 3,
	RGBA = 4
};

enum class InfoByte : uint8_t
{
	IsDct = 0x80,
	IsLong = 0x40,
};

struct SquashHeader
{
	uint32_t size_x;
	uint32_t size_y;

	ImageChannels channels;
};

struct CompressedBlock
{
	uint8_t infoByte;
	uint64_t table;

	int8_t data[BLOCK_SIZE * BLOCK_SIZE];
	uint8_t dataCount;
};

} // namespace sqh

#endif // INCLUDE_SQH_SQUASH_HEADER_HPP
