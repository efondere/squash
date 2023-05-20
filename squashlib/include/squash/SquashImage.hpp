/**
 * @file SquashImage.hpp
 * @author Eliot Fondere
 * @brief Main class for handling squash files (compression and decompression)
 *
 * @copyright Copyright (c) 2023 Eliot Fondere (MIT License)
 */

#ifndef INCLUDE_SQH_SQUASH_IMAGE_HPP
#define INCLUDE_SQH_SQUASH_IMAGE_HPP

#include "SquashHeader.hpp"
#include <string>

namespace sqh
{

class SquashImage
{
public:
	SquashImage(void* data, size_t dataSize, ImageChannels channels);
	explicit SquashImage(const std::string& file_path);

	~SquashImage();

	// prevent copying
	SquashImage(const SquashImage& other)            = delete;
	SquashImage& operator=(const SquashImage& other) = delete;

	bool open(const std::string& file_path);
	bool save(const std::string& file_path, bool overwrite=false);

	void* getData();

private:
	SquashHeader m_header;
	size_t       m_dataSize;
	void*        m_data;
};

} // sqh

#endif // INCLUDE_SQH_SQUASH_IMAGE_HPP
