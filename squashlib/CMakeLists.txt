add_library(squashlib
    include/squashlib/math/Matrix.hpp
    include/squashlib/math/math.hpp
    include/squashlib/squash.hpp
    include/squashlib/squash/squash_enc.hpp
    include/squashlib/squash/squash_dec.hpp
    include/squashlib/squash/SquashImage.hpp
    include/squashlib/squash/SquashHeader.hpp
    include/squashlib/utils/Buffer.hpp
    include/squashlib/utils/BufferD.hpp
    include/squashlib/math/MatrixD.hpp

    src/squashlib/squash/SquashImage.cpp
)

target_include_directories(squashlib
    PUBLIC
        thirdparty/stb/include
        include
)
