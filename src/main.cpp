/**
 * @file main.cpp
 * @author Eliot Fondere
 * @brief
 *
 * @copyright Copyright (c) 2023 Eliot Fondere (MIT License)
 */

#include "squash.hpp"

int main()
{
    float data[3] = {0.f, 2.f, 3.f};

	sqh::SquashImage image(&data, 3, sqh::ImageChannels::RGB);
	image.save("./image.sqh");
	sqh::SquashImage image2("./image.sqh");

	return 0;
}
