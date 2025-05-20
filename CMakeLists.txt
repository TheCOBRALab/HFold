cmake_minimum_required(VERSION 3.10)
project(HFold
    LANGUAGES C CXX
    DESCRIPTION "RNA folding prediction tool"
)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
option(ENABLE_TESTS "Build unit tests" ON)

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

# Default to Release build if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Get Raw user input for library path
set(_raw_user_input_lib "${VIENNARNA_LIB}") # user input using -D flag
if (NOT _raw_user_input_lib)
    set(_raw_user_input_lib "$ENV{VIENNARNA_LIB}") # user input using environment variable
endif()

# Set path to the cache variable VIENNARNA_LIB
if(_raw_user_input_lib)
    set(VIENNARNA_LIB
        "${_raw_user_input_lib}"
        CACHE PATH "Path to ViennaRNA headers" FORCE)
endif()

# Locate ViennaRNA library
find_library(VIENNARNA_LIBRARY
    NAMES vrna RNA
    HINTS ${VIENNARNA_LIB}
    PATHS
        /usr/local/lib/
        /usr/lib
        /opt/homebrew/lib
        /opt/vienna/lib
        /usr/local/lib64
        /usr/lib64
    DOC "Path to ViennaRNA library"
)

# ---------------- Clean include DIR path ----------------
# Get user input for include directory
set(_raw_user_input_dir "${VIENNARNA_INCLUDE_DIR}") # user input using -D flag
if (NOT _raw_user_input_dir)
    set(_raw_user_input_dir "$ENV{VIENNARNA_INCLUDE_DIR}") # user input using environment variable
endif()

set(_clean_user_input_dir "${_raw_user_input_dir}")

if(_clean_user_input_dir MATCHES "(.*/ViennaRNA)/pair_mat.h$")
    get_filename_component(_clean_user_input_dir "${_clean_user_input_dir}" DIRECTORY)  # strip pair_mat.h
    get_filename_component(_clean_user_input_dir "${_clean_user_input_dir}" DIRECTORY)  # strip ViennaRNA
elseif(_clean_user_input_dir MATCHES ".*/ViennaRNA/?$")
    get_filename_component(_clean_user_input_dir "${_clean_user_input_dir}" DIRECTORY) # strip ViennaRNA
endif()

# Set the cleaned path to the cache variable VIENNARNA_INCLUDE_DIR
if(_clean_user_input_dir)
    set(VIENNARNA_INCLUDE_DIR
        "${_clean_user_input_dir}"
        CACHE PATH "Path to ViennaRNA headers" FORCE)
endif()

# Search for the pair_mat.h file
find_path(VIENNARNA_INCLUDE_DIRECTORY
          NAMES  ViennaRNA/pair_mat.h
          HINTS "${VIENNARNA_INCLUDE_DIR}" # first try the cleaned path
          PATHS
                /usr/local/include
                /usr/include
                /opt/homebrew/include
                /opt/vienna/include
          DOC   "Path to ViennaRNA headers")

message(STATUS "VIENNARNA_LIBRARY = ${VIENNARNA_LIBRARY}")
message(STATUS "VIENNARNA_INCLUDE_DIRECTORY = ${VIENNARNA_INCLUDE_DIRECTORY}")



if(NOT VIENNARNA_LIBRARY OR NOT VIENNARNA_INCLUDE_DIRECTORY)
message(FATAL_ERROR "
❌ ViennaRNA library not found!

To fix this, either:

1. Set environment variables:
   export VIENNARNA_LIB=/path/to/libRNA.a
   export VIENNARNA_INCLUDE_DIR=/path/to/include

2. Or pass them to CMake directly:
   cmake -S . -B build \\
     -DVIENNARNA_LIB=/path/to/libRNA.a \\
     -DVIENNARNA_INCLUDE_DIR=/path/to/include

Tip: To find the files, try:
   find /usr/local /home -name libRNA.a
   find /usr/local /home -name pair_mat.h
")

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
target_include_directories(HFoldLib PRIVATE ${VIENNARNA_INCLUDE_DIRECTORY})
target_link_libraries(HFoldLib PUBLIC ${VIENNARNA_LIBRARY})
target_compile_options(HFoldLib PRIVATE ${COMMON_COMPILE_OPTIONS})

# Build the main executable
add_executable(HFold src/main.cpp)
target_sources(HFold PRIVATE ${hfold_sources}) # Re-add for Conda packaging (doesn't work without this for some reason)
target_include_directories(HFold PRIVATE ${VIENNARNA_INCLUDE_DIRECTORY}) 
target_link_libraries(HFold PRIVATE HFoldLib)
target_compile_options(HFold PRIVATE ${COMMON_COMPILE_OPTIONS})

# Install rule (For Conda Packaging)
install(TARGETS HFold RUNTIME DESTINATION bin)