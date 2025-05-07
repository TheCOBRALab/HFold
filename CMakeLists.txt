cmake_minimum_required(VERSION 3.9)
project(HFold
    LANGUAGES C CXX
    DESCRIPTION "RNA folding prediction tool"
)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/googletest/CMakeLists.txt")
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

# Default to Release build if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Locate ViennaRNA library
find_library(VIENNARNA_LIBRARY
    NAMES vrna RNA
    PATHS
        $ENV{VIENNARNA_LIB}
        /usr/local/lib
        /usr/lib
        /opt/homebrew/lib
        /opt/vienna/lib
        /usr/local/lib64
        /usr/lib64
    DOC "Path to ViennaRNA library"
)

message(STATUS "VIENNARNA_LIBRARY = ${VIENNARNA_LIBRARY}")

if(VIENNARNA_LIBRARY)
    message(STATUS "Found system ViennaRNA library at: ${VIENNARNA_LIBRARY}")
else()
    message(FATAL_ERROR "ViennaRNA library not found. Please install ViennaRNA and try again.")
endif()

# Core library sources (excluding main.cpp)
set(hfold_sources
    src/W_final.cpp 
    src/pseudo_loop.cpp 
    src/HFold.cpp
    src/cmdline.cpp
    src/Result.cpp	
    src/s_energy_matrix.cpp	
    src/Hotspot.cpp
    src/sparse_tree.cpp
)

# Common compile options
set(COMMON_COMPILE_OPTIONS
    $<$<CONFIG:Debug>:-g -Wall -Wextra> # cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
    $<$<CONFIG:Release>:-O3 -flto>      # cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
)

# Build the core logic as a library (Useful for testing)
add_library(HFoldLib ${hfold_sources})
target_link_libraries(HFoldLib PUBLIC ${VIENNARNA_LIBRARY})
target_compile_options(HFoldLib PRIVATE ${COMMON_COMPILE_OPTIONS})

# Build the main executable
add_executable(HFold src/main.cpp)
target_link_libraries(HFold PRIVATE HFoldLib)
target_compile_options(HFold PRIVATE ${COMMON_COMPILE_OPTIONS})

# Install rule (For Conda Packaging)
install(TARGETS HFold RUNTIME DESTINATION bin)