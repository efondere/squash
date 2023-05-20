/**
 * @file SquashImage.cpp
 * @author Eliot Fondere
 * @brief
 *
 * @copyright Copyright (c) 2023 Eliot Fondere (MIT License)
 */

#include "../../include/squash/SquashImage.hpp"
#include <fstream>

namespace sqh
{

SquashImage::SquashImage(void *data, size_t data_size, ImageChannels channels)
: m_header{128, 128, channels}
, m_dataSize(data_size)
, m_data(nullptr)
{
	m_data = new float[data_size];
	memcpy(m_data, data, data_size * sizeof(float));
}

SquashImage::SquashImage(const std::string &file_path)
: m_header{0, 0, ImageChannels::Grey}
, m_dataSize(0)
, m_data(nullptr)
{
	open(file_path);
}

bool SquashImage::open(const std::string &file_path)
{
	std::ifstream stream(file_path, std::ios::binary);
	stream.read(reinterpret_cast<char*>(&m_header), sizeof(m_header));
	stream.read(reinterpret_cast<char*>(&m_dataSize), sizeof(size_t));
	m_data = new float[m_dataSize];
	stream.read(reinterpret_cast<char*>(m_data), m_dataSize * sizeof(float));
	stream.close();

	return true;
}

bool SquashImage::save(const std::string &file_path, bool overwrite)
{
	std::ofstream stream(file_path, std::ios::binary);
	stream.write(reinterpret_cast<const char*>(&m_header), sizeof(m_header));
	stream.write(reinterpret_cast<const char*>(&m_dataSize), sizeof(size_t));
	stream.write(reinterpret_cast<const char*>(m_data), m_dataSize * sizeof(float));
	stream.close();

	return true;
}

SquashImage::~SquashImage()
{
	free(m_data);
}

void *SquashImage::getData()
{
	return m_data;
}

} // namespace sqh
