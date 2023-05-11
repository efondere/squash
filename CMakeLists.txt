cmake_minimum_required(VERSION 3.23)
project(squash)

set(CMAKE_CXX_STANDARD 17)

add_executable(squash src/main.cpp src/squash.hpp src/squash/squash_enc.hpp src/squash/squash_dec.hpp src/squash/SquashImage.cpp src/squash/SquashImage.hpp src/squash/SquashHeader.hpp src/math/Matrix.hpp src/math/math.hpp)

target_include_directories(squash PRIVATE ${PROJECT_SOURCE_DIR}/thirdparty/stb/include)
target_include_directories(squash PRIVATE ${PROJECT_SOURCE_DIR}/thirdparty/argparse/include)
