cmake_minimum_required(VERSION 3.9)

include(GoogleTest)

set(TEST_TARGET HFoldTests)

set(SOURCES
    HFoldTests.cpp
)

add_executable(${TEST_TARGET} ${SOURCES})
target_link_libraries(${TEST_TARGET} PUBLIC
    gtest_main
    HFoldLib
)

gtest_discover_tests(${TEST_TARGET}
    DISCOVERY_MODE PRE_TEST
    EXTRA_ARGS --gtest_color=yes
)

