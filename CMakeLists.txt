cmake_minimum_required(VERSION 3.9)
project(HFold
    LANGUAGES C CXX
    DESCRIPTION "RNA folding prediction tool"
)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

add_subdirectory(googletest)
add_subdirectory(test)

# Default to Release build if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# TODO: Find a more universal way to find the ViennaRNA library
# Locate ViennaRNA
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

find_path(VIENNARNA_INCLUDE_DIR
    NAMES ViennaRNA/fold.h
    PATHS
        $ENV{VIENNARNA_INC}
        /usr/local/include
        /usr/include
        /opt/homebrew/include
        /opt/vienna/include
        /usr/local/include/ViennaRNA
        /usr/include/ViennaRNA
    DOC "Path to ViennaRNA headers"
)

message(STATUS "VIENNARNA_LIBRARY = ${VIENNARNA_LIBRARY}")
message(STATUS "VIENNARNA_INCLUDE_DIR = ${VIENNARNA_INCLUDE_DIR}")

if(VIENNARNA_LIBRARY AND VIENNARNA_INCLUDE_DIR)
    message(STATUS "Found system ViennaRNA library at: ${VIENNARNA_LIBRARY}")
    message(STATUS "Found system ViennaRNA headers at: ${VIENNARNA_INCLUDE_DIR}")
    set(USE_SYSTEM_VIENNARNA TRUE)
else()
    message(FATAL_ERROR "ViennaRNA library and/or headers not found. Please install ViennaRNA and try again.")
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

# Build the core logic as a library
add_library(HFoldLib ${hfold_sources})
target_include_directories(HFoldLib PUBLIC ${VIENNARNA_INCLUDE_DIR})
target_link_libraries(HFoldLib PUBLIC ${VIENNARNA_LIBRARY})
target_compile_options(HFoldLib PRIVATE ${COMMON_COMPILE_OPTIONS})

# Build the main executable
add_executable(HFold src/main.cpp)
target_link_libraries(HFold PRIVATE HFoldLib)
target_compile_options(HFold PRIVATE ${COMMON_COMPILE_OPTIONS})

# Install rule (For Conda Packaging)
install(TARGETS HFold RUNTIME DESTINATION bin)
