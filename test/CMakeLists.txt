cmake_minimum_required(VERSION 3.9)

set(This ExampleTests)
set(SOURCES
    test_example.cpp
)

add_executable(${This} ${SOURCES})
target_link_libraries(${This} PUBLIC
    gtest_main
    HFoldLib
)

add_test(
    NAME ${This}
    COMMAND ${CMAKE_COMMAND} -E env
        "${CMAKE_CURRENT_BINARY_DIR}/${This}"
        --gtest_color=yes
        --gtest_filter=ExampleTests.*
)
