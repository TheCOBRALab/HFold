cmake_minimum_required(VERSION 3.9)

include(GoogleTest)

set(TEST_TARGET HFoldTests)

set(SOURCES
    tests/test_hfold.cpp
    tests/helper.cpp
)

add_executable(${TEST_TARGET} ${SOURCES})
target_link_libraries(${TEST_TARGET} PUBLIC
    gtest_main
    HFoldLib
)
target_include_directories(${TEST_TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/src)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/tests/input DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/tests/expected_output DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

gtest_discover_tests(${TEST_TARGET}
    DISCOVERY_MODE PRE_TEST
    EXTRA_ARGS --gtest_color=yes
)

