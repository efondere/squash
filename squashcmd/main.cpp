/**
 * @file main.cpp
 * @author Eliot Fondere
 * @brief
 *
 * @copyright Copyright (c) 2023 Eliot Fondere (MIT License)
 */

//https://github.com/nothings/stb/blob/master/stb_image_write.h
//https://github.com/nothings/stb/blob/master/stb_image.h
//https://github.com/p-ranav/argparse

#include "argparse/argparse.hpp"
#include "squashlib/squash/SquashImage.hpp"

#include <iostream>

int main(int argc, char* argv[])
{
	argparse::ArgumentParser program("squash", "0.1");

	program.add_argument("input")
		.required()
		.help("the input file to be compressed or decompressed");
	program.add_argument("-d")
		.implicit_value(true)
		.default_value(false)
		.help("decode the input file and save the decompressed data to the output file");
	program.add_argument("-c")
		.implicit_value(true)
		.default_value(false)
		.help("compress the input file and save the compressed data to the output file");
	program.add_argument("-o", "--output")
		.required()
		.metavar("output file")
		.help("the output file to save the compressed or decompressed data");

	try {
		program.parse_args(argc, argv);
	}
	catch (const std::runtime_error& err)
	{
		std::cerr << err.what() << std::endl;
		std::cerr << program;
		std::exit(1);
	}

	bool compress = program.get<bool>("-c");
	bool decompress = program.get<bool>("-d");

	if ((compress && decompress) || (!compress && !decompress))
	{
		std::cerr << "did you want to decompress (-d) or compress (-c)?" << std::endl;
		std::cerr << program;
		std::exit(1);
	}


	if (compress)
	{
		sqh::SquashImage::Quality = 0.8f;
		sqh::SquashImage img(program.get("input"));
		img.save(program.get("-o"), true);
	}
	else
	{
		sqh::SquashImage img(program.get("input"));
		img.save(program.get("-o"), true);
	}

	return 0;
}
