add_library(squashlib
    include/math/Matrix.hpp
    include/math/math.hpp
    include/squash.hpp
    include/squash/squash_enc.hpp
    include/squash/squash_dec.hpp
    include/squash/SquashImage.hpp
    include/squash/SquashHeader.hpp

    src/squash/SquashImage.cpp
)

target_include_directories(squashlib
    PUBLIC
        thirdparty/stb/include
        thirdparty/argparse/include
        include
)
