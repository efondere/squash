#include <squashlib/squash.hpp>

#include <filesystem>
#include <fstream>
#include <vector>

namespace fs = std::filesystem;

int main()
{
	// change this to use a different data set    \/
	auto root_path = fs::path("./data/TEXT_DATASET/");
	auto data_path = root_path / "data";
	auto sqh_out_path = root_path / "out";

	std::vector<uint64_t> compressed_sizes;
	std::vector<uint64_t> uncompressed_sizes;
	std::vector<float> file_average_error;

	sqh::SquashImage::Quality = 0.75f;

	for (const auto& entry : fs::directory_iterator(data_path))
	{
		auto sqh_file_path = sqh_out_path / (entry.path().filename().string() + ".sqh");

		sqh::SquashImage base_image(entry.path().string());
		base_image.save(sqh_file_path.string(), true);

		// extra info about file size, etc. is insignificant
		uncompressed_sizes.push_back(base_image.getHeader().size_x * base_image.getHeader().size_y * 3);
		compressed_sizes.push_back(fs::file_size(sqh_file_path));

		sqh::SquashImage compressed_image(sqh_file_path.string());

		double current_error = 0;

		for (size_t i = 0; i < base_image.getHeader().size_y; i++)
		{
			for (size_t j = 0; j < base_image.getHeader().size_x; j++)
			{
				for (size_t c = 0; c < 3; c++)
				{
					auto index = 3 * (base_image.getHeader().size_x * i + j) + c;
					auto dif = static_cast<float>(base_image.getData()[index]) - static_cast<float>(compressed_image.getData()[index]);

					current_error += dif * dif;
				}
			}
		}

		current_error /= base_image.getHeader().size_y * base_image.getHeader().size_x;
		std::cout << current_error << std::endl;

		file_average_error.push_back(static_cast<float>(current_error));
	}

	std::ofstream stats_file(root_path / "stats.txt", std::ios::out);

	for (auto cs : compressed_sizes)
	{
		stats_file << cs << " ";
	}
	stats_file << "\n";

	for (auto us : uncompressed_sizes)
	{
		stats_file << us << " ";
	}
	stats_file << "\n";

	for (auto er : file_average_error)
	{
		stats_file << er << " ";
	}
	stats_file << std::endl;

	stats_file.close();

	return 0;
}
