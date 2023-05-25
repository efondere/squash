add_executable(squashcmd
    main.cpp
)

target_link_libraries(squashcmd PUBLIC squashlib)

target_include_directories(squashcmd
    PRIVATE
        thirdparty/argparse/include
)
