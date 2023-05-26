add_library(squashlib
    include/squashlib/math/math.hpp
    include/squashlib/math/Matrix.hpp
    include/squashlib/squash/SquashHeader.hpp
    include/squashlib/squash/SquashImage.hpp
    include/squashlib/squash.hpp

    src/squashlib/squash/SquashImage.cpp
)

target_include_directories(squashlib
    PUBLIC
        thirdparty/stb/include
        include
)
