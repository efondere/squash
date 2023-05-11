/**
 * @file main.cpp
 * @author Eliot Fondere
 * @brief
 *
 * @copyright Copyright (c) 2023 Eliot Fondere (MIT License)
 */

#include "squash/squash_enc.hpp"
#include "squash/squash_dec.hpp"
#include "math/Matrix.hpp"
#include "math/math.hpp"

// https://github.com/nothings/stb/blob/master/stb_image_write.h
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>
#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>
//https://github.com/nothings/stb/blob/master/stb_image.h
//https://github.com/p-ranav/argparse

#include <argparse/argparse.hpp>

#include <iostream>
#include <cmath>

int main(int argc, char* argv[]) {
	argparse::ArgumentParser program("squash", "0.1");

	program.add_argument("eingabedatei")
		.required()
		.help("die Eingabedatei fur Kompression oder Dekodierung");
	program.add_argument("-d")
		.implicit_value(true)
		.default_value(false)
		.help("die Eingabedatei dekodieren und das entschlusseltes Bild auf dem Computer schreiben");
	program.add_argument("-c")
		.implicit_value(true)
		.default_value(false)
		.help("die Eingabedatei zukomprimieren und das komprimierte Bild auf dem Computer schreiben");
	program.add_argument("-o", "--output")
		.required()
		.metavar("AUSGABEDATEI")
		.help("wohin die Ausgabedatei geschreiben werden soll");

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
		std::cerr << "mochten Sie dekodiren (-d) oder zukomprimieren (-c)?" << std::endl;
		std::cerr << program;
		std::exit(1);
	}

	if (compress)
		sqh::encode(program.get("eingabedatei"), program.get("-o"));
	else
		sqh::decode(program.get("eingabedatei"), program.get("-o"));

	return 0;
}
