add_executable(HFoldTests
    tests/test_hfold.cc
    tests/helper.cc
)

target_link_libraries(HFoldTests PRIVATE
    HFoldLib
    GTest::gtest_main
)

include(GoogleTest)

# The tests use paths relative to the test binary directory. Copy the complete
# fixture directories so both input variants and expected output are present.
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/tests/input"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/tests/expected_output"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

gtest_discover_tests(HFoldTests
    DISCOVERY_TIMEOUT 30
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)
