/**
 * @file SquashImage.cpp
 * @author Eliot Fondere
 * @brief
 *
 * @copyright Copyright (c) 2023 Eliot Fondere (MIT License)
 */

#include <squashlib/squash/SquashImage.hpp>
#include <squashlib/math/math.hpp>

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image.h>
#include <stb/stb_image_write.h>

#include <filesystem>
#include <fstream>
#include <iostream>

namespace fs = std::filesystem;

namespace sqh
{

namespace
{

// --- DCT ---

// returns 1 if i == 0 and sqrt(2) otherwise
float delta_i(size_t i)
{
	if (i == 0)
		return 1.f;
	else
		return std::sqrt(2.f);
}

// gives the i, j entry of the C matrix for compression using DCT
float c_ij(size_t N, size_t i, size_t j)
{
	return (delta_i(i)/std::sqrt(static_cast<float>(N)))
	       * std::cos((static_cast<float>(i * (2 * j + 1)) * math::pi) / static_cast<float>(2 * N) );
}

// transformation matrix with DCT coefficients
const sqh::math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float> T_dct =
	math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float>::FromFunction([](size_t i, size_t j) {
		return c_ij(BLOCK_SIZE, i, j);
	});

// quantization table
float Q_dct_data_default[8][8] =
	{
		{10.f, 16.f, 22.f, 28.f, 34.f, 40.f, 46.f, 52.f},
		{16.f, 22.f, 28.f, 34.f, 40.f, 46.f, 52.f, 58.f},
		{22.f, 28.f, 34.f, 40.f, 46.f, 52.f, 58.f, 64.f},
		{28.f, 34.f, 40.f, 46.f, 52.f, 58.f, 64.f, 70.f},
		{34.f, 40.f, 46.f, 52.f, 58.f, 64.f, 70.f, 76.f},
		{40.f, 46.f, 52.f, 58.f, 64.f, 70.f, 76.f, 82.f},
		{46.f, 52.f, 58.f, 64.f, 70.f, 76.f, 82.f, 88.f},
		{52.f, 58.f, 64.f, 70.f, 76.f, 82.f, 88.f, 94.f}
	};

math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float> Q_dct_default = math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float>::FromArray(Q_dct_data_default);


// --- HAAR ---

float T_haar_data[BLOCK_SIZE][BLOCK_SIZE] =
	{
		{1.f/(2*std::sqrt(2.f)),  1.f/(2*std::sqrt(2.f)),  1.f/2,      0,  1/std::sqrt(2.f),             0,            0,            0},
		{1.f/(2*std::sqrt(2.f)),  1.f/(2*std::sqrt(2.f)),  1.f/2,      0, -1/std::sqrt(2.f),             0,            0,            0},
		{1.f/(2*std::sqrt(2.f)),  1.f/(2*std::sqrt(2.f)), -1.f/2,      0,             0,  1/std::sqrt(2.f),            0,            0},
		{1.f/(2*std::sqrt(2.f)),  1.f/(2*std::sqrt(2.f)), -1.f/2,      0,             0, -1/std::sqrt(2.f),            0,            0},
		{1.f/(2*std::sqrt(2.f)), -1.f/(2*std::sqrt(2.f)),      0,  1.f/2,             0,            0,  1/std::sqrt(2.f),            0},
		{1.f/(2*std::sqrt(2.f)), -1.f/(2*std::sqrt(2.f)),      0,  1.f/2,             0,            0, -1/std::sqrt(2.f),            0},
		{1.f/(2*std::sqrt(2.f)), -1.f/(2*std::sqrt(2.f)),      0, -1.f/2,             0,            0,            0,  1/std::sqrt(2.f)},
		{1.f/(2*std::sqrt(2.f)), -1.f/(2*std::sqrt(2.f)),      0, -1.f/2,             0,            0,            0, -1/std::sqrt(2.f)}
	};

float Q_haar_data_default[8][8] =
	{
		{ 8.f, 12.f, 16.f, 16.f, 24.f, 24.f, 24.f, 24.f},
		{12.f, 12.f, 16.f, 16.f, 24.f, 24.f, 24.f, 24.f},
		{16.f, 16.f, 24.f, 24.f, 32.f, 32.f, 32.f, 32.f},
		{16.f, 16.f, 24.f, 24.f, 32.f, 32.f, 32.f, 32.f},
		{24.f, 24.f, 32.f, 32.f, 38.f, 38.f, 38.f, 38.f},
		{24.f, 24.f, 32.f, 32.f, 38.f, 38.f, 38.f, 38.f},
		{24.f, 24.f, 32.f, 32.f, 38.f, 38.f, 38.f, 38.f},
		{24.f, 24.f, 32.f, 32.f, 38.f, 38.f, 38.f, 38.f}
	};



math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float> T_haar = math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float>::FromArray(T_haar_data).transpose();
math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float> Q_haar_default = math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float>::FromArray(Q_haar_data_default);


// --- UTILS ---
math::MatrixIndex zig_zag_indices[64] =
	{
		{0, 0},
		{0, 1}, {1, 0},
		{2, 0}, {1, 1}, {0, 2},
		{0, 3}, {1, 2}, {2, 1}, {3, 0},
		{4, 0}, {3, 1}, {2, 2}, {1, 3}, {0, 4},
		{0, 5}, {1, 4}, {2, 3}, {3, 2}, {4, 1}, {5, 0},
		{6, 0}, {5, 1}, {4, 2}, {3, 3}, {2, 4}, {1, 5}, {0, 6},
		{0, 7}, {1, 6}, {2, 5}, {3, 4}, {4, 3}, {5, 2}, {6, 1}, {7, 0},
		{7, 1}, {6, 2}, {5, 3}, {4, 4}, {3, 5}, {2, 6}, {1, 7},
		{2, 7}, {3, 6}, {4, 5}, {5, 4}, {6, 3}, {7, 2},
		{7, 3}, {6, 4}, {5, 5}, {4, 6}, {3, 7},
		{4, 7}, {5, 6}, {6, 5}, {7, 4},
		{7, 5}, {6, 6}, {5, 7},
		{6, 7}, {7, 6},
		{7, 7},
	};

math::MatrixIndex flatten_horiz(size_t index, size_t r, size_t c)
{
	return {index / c, index % c};
}

math::MatrixIndex zig_zag_flatten(size_t index)
{
	return zig_zag_indices[index];
}

constexpr uint64_t TABLE_BITMASK = uint64_t(1) << 63;
constexpr size_t DEFAULT_BLOCK_MEM_SIZE = BLOCK_SIZE * BLOCK_SIZE;
constexpr size_t OPTIMIZATION_ATTEMPTS = 3;
constexpr float LEARN_RATE = 0.25f;

} // namespace

double SquashImage::Quality = 0.5f;

SquashImage::SquashImage(std::string_view file_path)
: m_dctQTable(Q_dct_default)
, m_haarQTable(Q_haar_default)
{
	open(file_path);
}

SquashImage::~SquashImage()
{
	free();
}

bool SquashImage::open(std::string_view file_path)
{
	fs::path path(file_path);

	if (path.extension() == ".png" || path.extension() == ".PNG")
	{
		return open_png(file_path);
	}

	if (path.extension() == ".sqh")
	{
		return open_sqh(file_path);
	}

	std::cout << "[ERROR] (SquashImage): Unknown/unsupported file extension: \""
		<< path.extension() << "\"" << std::endl;
	return false;
}

bool SquashImage::save(std::string_view file_path, bool overwrite)
{
	fs::path path(file_path);

	if (path.extension() == ".png")
	{
		return save_png(file_path, overwrite);
	}

	if (path.extension() == ".sqh")
	{
		return save_sqh(file_path, overwrite);
	}

	std::cout << "[ERROR] (SquashImage): Unknown/unsupported file extension: \""
	          << path.extension() << "\"" << std::endl;
	return false;
}

bool SquashImage::open_png(std::string_view file_path)
{
	int width, height, channelCount;
	m_data = stbi_load(std::string(file_path).c_str(), &width, &height, &channelCount, 3);

	if (m_data == nullptr)
		return false;

	m_header.size_x = width;
	m_header.size_y = height;
	m_header.channels = ImageChannels::RGB; // always load & store RGB

	return true;
}

bool SquashImage::save_png(std::string_view file_path, bool overwrite)
{
	if (m_data == nullptr)
		return false;

	if (fs::exists(fs::path(file_path)) && !overwrite)
		return false;

	auto result = stbi_write_png(std::string(file_path).c_str(),
								 static_cast<int>(m_header.size_x), static_cast<int>(m_header.size_y),
								 3, m_data, 0);

	if (result == 0)
		return false;

	return true;
}

bool SquashImage::open_sqh(std::string_view file_path)
{
	std::ifstream input_file(std::string(file_path), std::ios::binary);

	uint32_t magic_number;
	input_file.read(reinterpret_cast<char*>(&magic_number), sizeof(uint32_t));

	if (magic_number != MAGIC_NUMBER)
	{
		std::cout << "[ERROR] (SquashImage): File \"" << file_path << "\" does not start with the correct magic number"
			<< std::endl;

		input_file.close();
		return false;
	}

	input_file.read(reinterpret_cast<char*>(&m_header), sizeof(SquashHeader));

	uint8_t q_data[BLOCK_SIZE][BLOCK_SIZE] = {};
	input_file.read(reinterpret_cast<char*>(q_data), sizeof(uint8_t) * BLOCK_SIZE * BLOCK_SIZE);
	m_dctQTable = math::Matrix<8, 8, uint8_t>::FromArray(q_data).asType<float>();
	input_file.read(reinterpret_cast<char*>(q_data), sizeof(uint8_t) * BLOCK_SIZE * BLOCK_SIZE);
	m_haarQTable = math::Matrix<8, 8, uint8_t>::FromArray(q_data).asType<float>();

	auto result = decompress(input_file);

	input_file.close();
	return result;
}

bool SquashImage::save_sqh(std::string_view file_path, bool overwrite)
{
	if  (m_data == nullptr)
		return false;

	std::ofstream output_file(std::string(file_path), std::ios::binary);

	output_file.write(reinterpret_cast<const char*>(&MAGIC_NUMBER), sizeof(uint32_t));
	output_file.write(reinterpret_cast<char*>(&m_header), sizeof(SquashHeader));

	auto result = compress(output_file);

	output_file.close();
	return result;
}

uint8_t* SquashImage::getData()
{
	return m_data;
}

const SquashHeader& SquashImage::getHeader()
{
	return m_header;
}

void SquashImage::free()
{
	::free(m_data);
}

math::Matrix<BLOCK_SIZE, BLOCK_SIZE, int8_t> SquashImage::transform_block(
	const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, uint8_t> &block,
    const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float> &transform_matrix,
    const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float> &quantization_matrix,
	math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float>* coefficients)
{
	// AKA f_tilde
	auto shifted_data = block.asType<float>([](uint8_t input) -> float {
		return static_cast<float>(input) - 128.f;
	});

	// AKA alpha_tilde (equation 12.9)
	auto transformed_data = transform_matrix.product(shifted_data.product(transform_matrix.transpose()));
	if (coefficients != nullptr)
		*coefficients = transformed_data;

	// AKA l
	auto quantized_data = (transformed_data / quantization_matrix) + 0.5f;
	return quantized_data.asType<int8_t>([](float input) {
		return static_cast<int8_t>(std::floorf(input));
	});
}

math::Matrix<8, 8, uint8_t> SquashImage::inverse_transform_block(
	std::ifstream &input_file, uint8_t info_byte,
    const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float> &transform_matrix,
    const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float> &quantization_matrix)
{
	auto block = decompress_block(input_file, info_byte);

	auto beta = (block.asType<float>() * quantization_matrix);

	auto f_bar = transform_matrix.transpose().product(beta.product(transform_matrix)) + 128.f;
	return f_bar.asType<uint8_t>([](float input) {
		return static_cast<uint8_t>(std::min(255.f, std::max(0.f, floorf(input))));
	});
}

math::Matrix<8, 8, uint8_t> SquashImage::test_inverse_transform_block(
	const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, int8_t>& input_block,
	const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float>& transform_matrix,
	const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float>& quantization_matrix)
{
	auto beta = (input_block.asType<float>() * quantization_matrix);

	auto f_bar = transform_matrix.transpose().product(beta.product(transform_matrix)) + 128.f;
	return f_bar.asType<uint8_t>([](float input) {
		return static_cast<uint8_t>(std::min(255.f, std::max(0.f, floorf(input))));
	});
}

void SquashImage::compress_block(
	const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, int8_t>& block_data, CompressedBlock& compressed_block,
	bool isHaar)
{
	compressed_block.infoByte = isHaar ? 0 : static_cast<uint8_t>(InfoByte::IsDct);

	// STEP 1: try short with zig-zag
	auto zig_zag_data = block_data.flatten(zig_zag_flatten);
	size_t endZerosCount = 0;
	int index = zig_zag_data.size() - 1;
	while (index >= 0 && zig_zag_data[index] == 0)
	{
		endZerosCount++;
		index--;
	}

	size_t extraZeros = 0;
	while (index >= 0)
	{
		if (zig_zag_data[index] == 0)
		{
			extraZeros++;
		}
		index--;
	}

	if (extraZeros <= 8)
	{
		// use short representation
		auto dataCount = static_cast<uint8_t>(64 - endZerosCount); // TODO: careful here! If the block is completely full, this will lead to problems...
		compressed_block.infoByte |= dataCount;
		compressed_block.table = 0;
		compressed_block.dataCount = dataCount;

		for (uint8_t i = 0; i < dataCount; i++)
		{
			compressed_block.data[i] = zig_zag_data[i];
		}

		return;
	}

	compressed_block.infoByte |= static_cast<uint8_t>(InfoByte::IsLong);

	// we could store this in zig-zag to avoid recomputing but whatever
	auto horizontal_data = block_data.flatten(flatten_horiz);

	for (auto value : horizontal_data)
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

math::Matrix<8, 8, int8_t> SquashImage::decompress_block(
	std::ifstream &input_file, uint8_t infoByte)
{
	math::Matrix<8, 8, int8_t> m;

	if (infoByte & static_cast<uint8_t>(InfoByte::IsLong))
	{
		// long block
		uint64_t table;
		input_file.read(reinterpret_cast<char*>(&table), sizeof(table));

		for (int k = 0; k < 64; k++)
		{

			auto index = std::div(k, 8);

			if (!(table & TABLE_BITMASK))
			{
				m.data[index.quot][index.rem] = 0;
			}
			else
			{
				int8_t value;
				input_file.read(reinterpret_cast<char*>(&value), sizeof(value));
				m.data[index.quot][index.rem] = value;
			}

			table <<= 1;
		}

		return m;
	}

	// short block
	uint8_t data_count = infoByte & 0x3F; // 63 (so all ones except the first two)
	for (int i = 0; i < data_count; i++)
	{
		int8_t value;
		input_file.read(reinterpret_cast<char*>(&value), sizeof(value));

		auto index = zig_zag_flatten(i);
		m.data[index.i][index.j] = value;
	}

	return m;
}

bool SquashImage::decompress(std::ifstream &input_file)
{
	auto x_div = std::div(static_cast<int64_t>(m_header.size_x), BLOCK_SIZE);
	auto y_div = std::div(static_cast<int64_t>(m_header.size_y), BLOCK_SIZE);
	uint32_t x_blocks = x_div.quot + (x_div.rem == 0 ? 0 : 1);
	uint32_t y_blocks = y_div.quot + (y_div.rem == 0 ? 0 : 1);

	free();
	m_data = reinterpret_cast<uint8_t*>(malloc(m_header.size_x * m_header.size_y * 3));

	for (int i = 0; i < y_blocks; i++) {
		for (int j = 0; j < x_blocks; j++) {
			for (int c = 0; c < 3; c++) {
				uint8_t info_byte = 0;
				input_file.read(reinterpret_cast<char*>(&info_byte), sizeof(info_byte));

				math::Matrix<BLOCK_SIZE, BLOCK_SIZE, uint8_t> f_bar;

				if (info_byte & static_cast<uint8_t>(InfoByte::IsDct))
				{
					// DCT
					f_bar = inverse_transform_block(input_file, info_byte, T_dct, m_dctQTable);
				}
				else
				{
					// HAAR
					f_bar = inverse_transform_block(input_file, info_byte, T_haar, m_haarQTable);
				}

				for (int k = 0; k < 8; k++) { // row
					for (int l = 0; l < 8; l++) { // col
						if ((8 * i + k >= m_header.size_y) || (8 * j + l >= m_header.size_x)) continue;
						m_data[3 * ((m_header.size_x * ((BLOCK_SIZE * i) + k)) + (BLOCK_SIZE * j) + l) + c] = f_bar.data[k][l];
					}
				}
			}
		}
	}

	return true;
}

bool SquashImage::compress(std::ofstream &output_file)
{
	//findOptimalQTables();

	auto Q_haar_data = m_haarQTable.asType<uint8_t>().flatten(flatten_horiz);
	auto Q_dct_data = m_dctQTable.asType<uint8_t>().flatten(flatten_horiz);
	output_file.write(reinterpret_cast<const char*>(Q_dct_data.data()), BLOCK_SIZE * BLOCK_SIZE);
	output_file.write(reinterpret_cast<const char*>(Q_haar_data.data()), BLOCK_SIZE * BLOCK_SIZE);

	auto x_div = std::div(static_cast<int64_t>(m_header.size_x), BLOCK_SIZE);
	auto y_div = std::div(static_cast<int64_t>(m_header.size_y), BLOCK_SIZE);
	uint32_t x_blocks = x_div.quot + (x_div.rem == 0 ? 0 : 1);
	uint32_t y_blocks = y_div.quot + (y_div.rem == 0 ? 0 : 1);

	// STATS
	double averageCompressionQuality = 0.0;

	for (int i = 0; i < y_blocks; i++) {
		for (int j = 0; j < x_blocks; j++) {
			for (int c = 0; c < 3; c++) {
				auto block = math::Matrix<BLOCK_SIZE, BLOCK_SIZE, uint8_t>::FromFunction(
					[this, i, j, c](size_t k, size_t l) -> uint8_t {
						if (8 * j + l >= m_header.size_x) return 128;
						if (8 * i + k >= m_header.size_y) return 128;

						return m_data[3 * ((m_header.size_x * ((BLOCK_SIZE * i) + k)) + (BLOCK_SIZE * j) + l) + c];
					});

				auto block_dct = transform_block(block, T_dct, m_dctQTable);
				auto block_haar = transform_block(block, T_haar, m_haarQTable);

				CompressedBlock compressed_dct{};
				CompressedBlock compressed_haar{};
				CompressedBlock* best_block;

				compress_block(block_dct, compressed_dct, false);
				auto dct_quality = computeCompressionQuality(
					block, test_inverse_transform_block(block_dct, T_dct, m_dctQTable),
					getCompressedSize(compressed_dct));
				compress_block(block_haar, compressed_haar, true);
				auto haar_quality = computeCompressionQuality(
					block, test_inverse_transform_block(block_haar, T_haar, m_haarQTable),
					getCompressedSize(compressed_haar));

				if (abs(Quality - haar_quality) < abs(Quality - dct_quality))
				{
					best_block = &compressed_haar;
					averageCompressionQuality += haar_quality;
				}
				else
				{
					best_block = &compressed_dct;
					averageCompressionQuality += dct_quality;
				}

				output_file.write(reinterpret_cast<const char*>(
									&best_block->infoByte), sizeof(best_block->infoByte));
				if (best_block->infoByte & static_cast<uint8_t>(InfoByte::IsLong))
				{
					output_file.write(reinterpret_cast<const char *>(
						                  &best_block->table), sizeof(best_block->table));
				}
				output_file.write(reinterpret_cast<const char *>(
					                  &best_block->data), sizeof(int8_t) * best_block->dataCount);
			}
		}
	}

	averageCompressionQuality /= static_cast<float>(y_blocks) * static_cast<float>(x_blocks) * 3.f;
	std::cout << "Compression success: " << averageCompressionQuality << " compared to requested: " << Quality << std::endl;

	return true;
}

size_t SquashImage::getCompressedSize(CompressedBlock& compressed_block)
{
	size_t totalSize = 1; // info byte

	if (compressed_block.infoByte & static_cast<uint8_t>(InfoByte::IsLong))
	{
		totalSize += 8;
	}

	totalSize += compressed_block.dataCount;

	return totalSize;
}

double SquashImage::computeCompressionQuality(const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, uint8_t> &original,
                                             const math::Matrix<BLOCK_SIZE, BLOCK_SIZE, uint8_t> &compressed,
                                             size_t compressed_size)
{

	// 1 is identical, 0 is "completely different"
	//auto block_resemblance = original.asType<double>().dot(compressed.asType<double>()) / 100000 / 50; // (0.25 - 0.90)
	auto block_resemblance = exp(-(original.asType<float>() - compressed.asType<float>()).norm() / 64.f);

	// 0 is impossible >= 1: means no / bad compression
	auto compression_ratio = static_cast<double>(compressed_size) / static_cast<double>(DEFAULT_BLOCK_MEM_SIZE);

	return sqrt(block_resemblance * block_resemblance + compression_ratio * compression_ratio);
}

void SquashImage::findOptimalQTables()
{

	// TO COMPUTE DIVERGENCE: Quality - compressionQuality
	// +ve: block is very nice and not very compressed
	// -ve: block is ugly and too compressed
	auto x_div = std::div(static_cast<int64_t>(m_header.size_x), BLOCK_SIZE);
	auto y_div = std::div(static_cast<int64_t>(m_header.size_y), BLOCK_SIZE);
	uint32_t x_blocks = x_div.quot + (x_div.rem == 0 ? 0 : 1);
	uint32_t y_blocks = y_div.quot + (y_div.rem == 0 ? 0 : 1);

	for (size_t attempt = 0; attempt < OPTIMIZATION_ATTEMPTS; attempt++)
	{
		math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float> averageHaarTransform;
		math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float> averageDctTransform;
		double averageDctQuality = 0.0;
		double averageHaarQuality = 0.0;

		for (int i = 0; i < y_blocks; i++) {
			for (int j = 0; j < x_blocks; j++) {
				for (int c = 0; c < 3; c++) {

					auto block = math::Matrix<BLOCK_SIZE, BLOCK_SIZE, uint8_t>::FromFunction(
						[this, i, j, c](size_t k, size_t l) -> uint8_t {
							if (8 * j + l >= m_header.size_x) return 128;
							if (8 * i + k >= m_header.size_y) return 128;

							return m_data[3 * ((m_header.size_x * ((BLOCK_SIZE * i) + k)) + (BLOCK_SIZE * j) + l) + c];
						});

					math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float> transformed_block_haar;
					auto block_haar = transform_block(block, T_haar, m_haarQTable, &transformed_block_haar);
					math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float> transformed_block_dct;
					auto block_dct = transform_block(block, T_dct, m_dctQTable, &transformed_block_dct);

					CompressedBlock compressed_haar{};
					CompressedBlock compressed_dct{};

					compress_block(block_haar, compressed_haar, true);
					auto haar_quality = computeCompressionQuality(
						block, test_inverse_transform_block(block_haar, T_haar, m_haarQTable),
						getCompressedSize(compressed_haar));
					compress_block(block_dct, compressed_dct, true);
					auto dct_quality = computeCompressionQuality(
						block, test_inverse_transform_block(block_dct, T_dct, m_dctQTable),
						getCompressedSize(compressed_dct));

					if (abs(Quality - haar_quality) <= abs(Quality - dct_quality))
					{
						// use haar
						averageHaarTransform = averageHaarTransform + transformed_block_haar;
						averageHaarQuality += haar_quality;
					}
					else
					{
						// use dct
						averageDctTransform = averageDctTransform + transformed_block_dct;
						averageDctQuality += dct_quality;
					}
				}
			}
		}

		auto total_blocks = static_cast<float>(y_blocks * x_blocks * 3);
		averageDctTransform = averageDctTransform / total_blocks;
		averageHaarTransform = averageHaarTransform / total_blocks;
		averageDctQuality /= total_blocks;
		averageHaarQuality /= total_blocks;
		std::cout << averageHaarQuality << std::endl;

		// process haar:
		if ((Quality - averageHaarQuality) < 0.0)
		{
			std::cout << "HIGHER HAAR" << std::endl;

			//m_haarQTable = m_haarQTable * (1.f + LEARN_RATE);

			//m_haarQTable = (m_haarQTable + averageHaarTransform.asType2<float>([&attempt](float value)
            //   {
            //       if (abs(value) <= 0.000001f) return 50.f;
            //       else return fmin(50.f, (LEARN_RATE / abs(value) / static_cast<float>(attempt)));
            //   })).asType<float>([](float value){return fmax(1.f, value);});
		}
		else if ((Quality - averageHaarQuality) > 0.0)
		{

			//m_haarQTable = m_haarQTable / (1.f + LEARN_RATE);

			//m_haarQTable = (m_haarQTable - averageHaarTransform.asType2<float>([&attempt](float value)
	        //   {
		    //       if (abs(value) <= 0.000001f) return 50.f;
		    //       else return fmin(50.f, (LEARN_RATE / abs(value) / static_cast<float>(attempt)));
	        //   })).asType<float>([](float value){return fmax(1.f, value);});
		}

		// process dct:
		if ((Quality - averageDctQuality) < 0.0)
		{
			std::cout << "HIGHER DCT" << std::endl;

			//m_dctQTable = (m_dctQTable + averageDctTransform.asType2<float>([&attempt](float value)
            //   {
            //       if (abs(value) <= 0.000001f) return 50.f;
            //       else return fmin(50.f, (LEARN_RATE / abs(value) / static_cast<float>(attempt)));
            //   })).asType<float>([](float value){return fmax(1.f, value);});

			m_dctQTable = m_dctQTable * (1.f + LEARN_RATE);
		}
		else if ((Quality - averageDctQuality) > 0.0)
		{
			std::cout << "LOWER DCT" << std::endl;

			m_dctQTable = m_dctQTable / (1.f + LEARN_RATE);

			//m_dctQTable = (m_dctQTable - averageDctTransform.asType2<float>([&attempt](float value)
            //   {
            //       if (abs(value) <= 0.000001f) return 10.f;
            //       else return fmin(10.f, (LEARN_RATE / abs(value) / static_cast<float>(attempt)));
            //   })).asType<float>([](float value){return fmax(8.f, value);});
		}
	}

	std::ofstream output_file("./output_matrices.txt", std::ios::out);
	m_haarQTable.asType<int>().toStream(output_file);
	m_dctQTable.asType<int>().toStream(output_file);
	output_file.close();
}

/*
void SquashImage::findOptimalDCTQTable()
{
	auto x_div = std::div(static_cast<int64_t>(m_header.size_x), BLOCK_SIZE);
	auto y_div = std::div(static_cast<int64_t>(m_header.size_y), BLOCK_SIZE);
	uint32_t x_blocks = x_div.quot + (x_div.rem == 0 ? 0 : 1);
	uint32_t y_blocks = y_div.quot + (y_div.rem == 0 ? 0 : 1);

	size_t attempts = 0;
	while(attempts < OPTIMIZATION_ATTEMPTS)
	{
		// STATS
		double averageCompressionQuality = 0.0;

		for (int i = 0; i < y_blocks; i++) {
			if (attempts >= OPTIMIZATION_ATTEMPTS)
				break;

			for (int j = 0; j < x_blocks; j++) {
				if (attempts >= OPTIMIZATION_ATTEMPTS)
					break;

				for (int c = 0; c < 3; c++) {
					if (attempts >= OPTIMIZATION_ATTEMPTS)
						break;

					attempts ++;

					auto block = math::Matrix<BLOCK_SIZE, BLOCK_SIZE, uint8_t>::FromFunction(
						[this, i, j, c](size_t k, size_t l) -> uint8_t {
							if (8 * j + l >= m_header.size_x) return 128;
							if (8 * i + k >= m_header.size_y) return 128;

							return m_data[3 * ((m_header.size_x * ((BLOCK_SIZE * i) + k)) + (BLOCK_SIZE * j) + l) + c];
						});

					math::Matrix<BLOCK_SIZE, BLOCK_SIZE, float> transformed_block;
					auto block_dct = transform_block(block, T_dct, m_dctQTable, &transformed_block);

					CompressedBlock compressed_dct{};

					compress_block(block_dct, compressed_dct, true);
					auto dct_quality = computeCompressionQuality(
						block, test_inverse_transform_block(block_dct, T_dct, m_dctQTable),
						getCompressedSize(compressed_dct));

					if ((Quality - dct_quality) < 0.0)
					{
						m_dctQTable = (m_dctQTable + transformed_block.asType<float>([](float value)
						                                                               {
							                                                               if (abs(value) <= 0.00001f) return 1.f;
							                                                               else return fmin(1.f, LEARN_RATE / abs(value));
						                                                               })).asType<float>([](float value){return fmax(1.f, value);});
					}
					else if ((Quality - dct_quality) > 0.0)
					{
						m_dctQTable = (m_dctQTable - transformed_block.asType<float>([](float value)
						                                                               {
							                                                               if (abs(value) <= 0.00001f) return 1.f;
							                                                               else return fmin(1.f, LEARN_RATE / abs(value));
						                                                               })).asType<float>([](float value){return fmax(1.f, value);});
					}
				}
			}
		}
	}
}*/

} // namespace sqh
