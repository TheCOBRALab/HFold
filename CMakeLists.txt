cmake_minimum_required(VERSION 3.15)
project(HFold)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_TESTS "Build unit tests" ON)
option(BUILD_PYTHON "Build Python bindings using pybind11" OFF)

# Directory for parameters
set(PARAMS_DIR "${CMAKE_SOURCE_DIR}/params")  # or CMAKE_INSTALL_PREFIX based path

if (NOT ENABLE_TESTS)
    message(STATUS "Unit tests are disabled. To enable them, set -DBUILD_TESTING=ON")
elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/googletest/CMakeLists.txt")
    enable_testing()
    add_subdirectory(test)
    add_subdirectory(test/googletest)
else()
    message(WARNING "
============================================================
 GoogleTest submodule not found. Unit tests will be disabled.
------------------------------------------------------------
To enable unit testing, please do one of the following:

  1. If you already cloned the repository:
     Run the following command to initialize submodules:
         git submodule update --init --recursive

  2. If you haven't cloned the repository yet:
     Clone it with submodules using:
         git clone --recursive <repo-url>

After initializing the submodule, re-run CMake to enable tests.
============================================================")

endif()

# Default build type if none is specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# -------------------------------------------------------------------
#                             C++ Flags
# -------------------------------------------------------------------

# Debug configuration
if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "Strict")
    message(STATUS "Building in Debug mode")
    add_compile_options(
        -g
        -O0
        -Wall
        -Wextra
        # -Wpedantic
        # -Wshadow
        # -Wsign-compare
        # -Wconversion
        # -Wsign-conversion
        -fsanitize=address,undefined
        -fno-sanitize-recover=all
        -fno-omit-frame-pointer
    )

    if(CMAKE_BUILD_TYPE STREQUAL "Strict")
        add_compile_options(-Werror)
    endif()

    add_link_options(
        -fsanitize=address,undefined
        -fno-sanitize-recover=all
    )
else()
    # Normal (Release or other non-Debug builds)
    message(STATUS "Building in Release/Normal mode")
    add_compile_options(
        -O3
        -march=native
        -Wall
        -Wextra
    )
endif()

set(SOURCE
  src/main.cc
  src/pseudo_loop_trace.cc
  src/pseudo_loop.cc 
  src/HFold.cc
  src/cmdline.cc
  src/Result.cc	
  src/s_energy_matrix.cc	
  src/Hotspot.cc
  src/sparse_tree.cc
  src/SHAPE.cc
)

set(vienna_SOURCE
    src/ViennaRNA/params/default.cc
    src/ViennaRNA/params/params.cc
    src/ViennaRNA/params/io.cc
    src/ViennaRNA/model.cc 
    src/ViennaRNA/utils.cc
)

add_library(RNA STATIC ${vienna_SOURCE})

target_include_directories(RNA PRIVATE .)

add_executable(HFold ${SOURCE})
target_sources(HFold PRIVATE ${SOURCE}) # Re-add for Conda packaging (doesn't work without this for some reason)
# Pass it as a compile definition
target_compile_definitions(HFold PRIVATE 
    PARAMS_DIR="${PARAMS_DIR}"
)

# Install CParty into the bin directory
install(TARGETS HFold RUNTIME DESTINATION bin)

target_link_libraries(HFold PRIVATE RNA)

include_directories(src)

if (BUILD_PYTHON)
    add_subdirectory(pylib)
endif()
