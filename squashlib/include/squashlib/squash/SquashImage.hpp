/**
 * @file SquashImage.hpp
 * @author Eliot Fondere
 * @brief Main class for handling squash files (compression and decompression)
 *
 * @copyright Copyright (c) 2023 Eliot Fondere (MIT License)
 */

#ifndef INCLUDE_SQH_SQUASH_IMAGE_HPP
#define INCLUDE_SQH_SQUASH_IMAGE_HPP

#include <squashlib/squash/SquashHeader.hpp>
#include <squashlib/math/Matrix.hpp>
#include <array>
#include <string>

namespace sqh
{

class SquashImage
{
public:
	static double Quality;

	explicit SquashImage(std::string_view file_path);
	~SquashImage();

	// prevent copying
	SquashImage(const SquashImage& other)            = delete;
	SquashImage& operator=(const SquashImage& other) = delete;

	bool open(std::string_view file_path);
	bool save(std::string_view file_path, bool overwrite);

	bool open_png(std::string_view file_path);
	bool save_png(std::string_view file_path, bool overwrite=false);

	bool open_sqh(std::string_view file_path);
	bool save_sqh(std::string_view file_path, bool overwrite=false);

	uint8_t* getData();
	const SquashHeader& getHeader();

	const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float>& getDCTQTable() const;
	const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float>& getHaarQTable() const;

	void free();

private:
	static math::Matrix<BLOCK_SIZE, BLOCK_SIZE, int8_t> transform_block(
		const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, uint8_t>& block,
		const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float>& transform_matrix,
		const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float>& quantization_matrix,
		math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float>* coefficients = nullptr);

	static math::Matrix<8, 8, uint8_t> inverse_transform_block(
		std::ifstream& input_file, uint8_t info_byte,
		const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float>& transform_matrix,
		const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float>& quantization_matrix);

	static math::Matrix<8, 8, uint8_t> test_inverse_transform_block(
		const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, int8_t>& input_block,
		const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float>& transform_matrix,
		const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float>& quantization_matrix);

	static void compress_block(
		const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, int8_t>& block_data, CompressedBlock& compressed_block,
		bool isHaar);
	static math::Matrix<8, 8, int8_t> decompress_block(
		std::ifstream& input_file, uint8_t infoByte);

	bool decompress(std::ifstream& input_file);
	bool compress(std::ofstream& output_file);

	static size_t getCompressedSize(CompressedBlock& compressed_block);

	static double computeCompressionQuality(
		const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, uint8_t>& original,
		const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, uint8_t>& compressed,
		size_t compressedSize);

	void findOptimalQTables();

	SquashHeader m_header{};
	uint8_t*     m_data = nullptr;

	math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float> m_dctQTable;
	math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float> m_haarQTable;
};

} // sqh

#endif // INCLUDE_SQH_SQUASH_IMAGE_HPP
