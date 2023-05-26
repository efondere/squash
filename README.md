# squash
Compressed image format

Recommended extension: .sqh

## Building the project

This project uses CMake as a build system. Therefore, it is recommended to open it in an editor supporting CMake, such
as [CLion](https://www.jetbrains.com/clion/) (free educational licenses) or
[Visual Studio](https://learn.microsoft.com/en-us/cpp/build/cmake-projects-in-visual-studio?view=msvc-170). One may have
to install CMake on their computer [from here](https://cmake.org/). This project can also be generated (that is, CMake
will generate project files for a specific IDE / Compiler) manually by opening the CMake GUI. This application comes with
CMake when installing it.

## Usage

This project builds a library and two executables:
- _squashlib_: this is the library used by the two executables. It contains the compression algorithms
- _squashcmd_: this is a simple command-line tool which allows to compress or decompress individual images. For more
info, run `squashcmd.exe --help` in a terminal.
- _squashtest_: this executable will go through all files in a directory to and compress them to test the efficiency of
the compression algorithm. See bellow for usage.

## Squashtest

Test images curated by us can be found in the
[releases](https://github.com/espressjo/ommtools/releases) page of the repository.

This executable will loop through all .png files of a given directory, compress them, save the compressed images in
squash files, and generate a text file containing statistics about the compressed images. Statistics collected are:
- compressed file size
- hypothetical fully uncompressed data size (ie. width * height * number of color channels)
- the norm of the difference between the matrices representing the uncompressed and compressed images, divided by the
size of the image.

The folder structure for test data is the following:
```
└─ ROOT
    ├─ data
    ├─ out 
    └─ stats.txt
```
Where `ROOT` will usually take the name of the dataset. The data folder contains all the images and the out folder will
be filled with all squash files resulting from the compression. Finally, the stats.txt files will contain the statistics
collected during the test. Each line contains the values for each collected statistic, where the data for each image is
separated by a space. The path to the root folder must be changed on line 12 of the `main.cpp` in the squashtest
directory of this project. Note that the path is relative to where the executable (squashtest.exe) is run. This usually
depends on the IDE, but can be determined once the project has been built completely at least once.
