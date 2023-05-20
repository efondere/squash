cmake_minimum_required(VERSION 3.23)
project(squash CXX)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(squashlib)
add_subdirectory(squashcmd)
add_subdirectory(squashtest)
