#include <filesystem>

#include <squashlib/squash.hpp>
#include <squashlib/math/MatrixD.hpp>

namespace fs = std::filesystem;

constexpr float learn_rate = 1.f;

int main()
{
	std::string path = "./data/DIV2K/";
	double prev_error = -1.0;

	//sqh::priv::q = sqh::math::Matrix<8, 8, float>::FromFunction([] (size_t, size_t){ return 20.f; });

	for (const auto& entry : fs::directory_iterator(path))
	{
		auto sqh_path = std::string("./data/DIV2K/out_sqh/") + entry.path().filename().string() + ".sqh";
		auto output_path = std::string("./data/DIV2K/out/") + entry.path().filename().string();
		auto avDCT = sqh::encode2(entry.path().string(), sqh_path);
		sqh::decode(sqh_path, output_path);

		double compression_ratio = static_cast<double>(fs::file_size(fs::path(sqh_path))) / static_cast<double>(entry.file_size());

		int w1, h1, n1, w2, h2, n2;
		uint8_t* data1 = stbi_load(entry.path().string().c_str(), &w1, &h1,
		                                &n1, 1);
		uint8_t* data2 = stbi_load(output_path.c_str(), &w2, &h2,
		                                &n2, 1);

		sqh::math::MatrixD<uint8_t> m1(w1, h1);
		m1.fromArray(data1);
		sqh::math::MatrixD<uint8_t> m2(w1, h1);
		m2.fromArray(data2);

		auto dif = m1 - m2;
		double norm = 1000 * dif.norm() / (w1 * h1);

		stbi_image_free(data1);
		stbi_image_free(data2);


		auto error = sqrt(compression_ratio * compression_ratio + norm * norm);
		std::cout << entry.path().filename().string() << compression_ratio << " + " << norm << " = " << error << std::endl;

		if (prev_error >= 0)
		{
			if (error > prev_error) {
				sqh::priv::q = sqh::priv::q + (-learn_rate / avDCT);
			}
			else
				sqh::priv::q = sqh::priv::q + (learn_rate / avDCT);
		}

		prev_error = error;

		// 1) process first image
		// 2) return the average DCT coefficients
		// 3) calculate norm of angle and compression ratio vector
		// 3) add learn_speed/avDCT to q and pass new q for second image
		// 4) calculate the new norm that we are trying to minimize
		// 5) decrease? -> add calculate new avDCT and add learn_speed/avDCT to quantization table
		// 6) increase? -> subtract learn_speed/2 * 1/avDCT
	}

	sqh::priv::q.print();

	return 0;
}
