cmake_minimum_required(VERSION 3.16)

project(HFold VERSION 1.16.0 LANGUAGES C CXX)

include(GNUInstallDirs)

option(ENABLE_TESTS "Build unit tests" ON)
option(BUILD_PYTHON "Build Python bindings using pybind11" OFF)
option(HFOLD_ENABLE_SANITIZERS
       "Enable AddressSanitizer and UndefinedBehaviorSanitizer in Debug builds"
       ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                 Debug Release RelWithDebInfo MinSizeRel)
endif()

# Shared compiler and linker settings. Because this is an INTERFACE target,
# executables that use HFoldLib also receive the sanitizer link flags required
# by the instrumented static libraries.
add_library(HFoldBuildOptions INTERFACE)
target_compile_features(HFoldBuildOptions INTERFACE cxx_std_17)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(HFoldBuildOptions INTERFACE
        -Wall
        -Wextra
        "$<$<CONFIG:Debug>:-O0>"
    )

    if(HFOLD_ENABLE_SANITIZERS)
        target_compile_options(HFoldBuildOptions INTERFACE
            "$<$<CONFIG:Debug>:-fsanitize=address,undefined;-fno-sanitize-recover=all;-fno-omit-frame-pointer>"
        )
        target_link_options(HFoldBuildOptions INTERFACE
            "$<$<CONFIG:Debug>:-fsanitize=address,undefined;-fno-sanitize-recover=all>"
        )
    endif()
endif()

# Embedded ViennaRNA subset used by HFold.
add_library(RNA STATIC
    src/ViennaRNA/params/default.cc
    src/ViennaRNA/params/params.cc
    src/ViennaRNA/params/io.cc
    src/ViennaRNA/model.cc
    src/ViennaRNA/utils.cc
)

target_include_directories(RNA PUBLIC
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/src"
)

target_link_libraries(RNA PUBLIC HFoldBuildOptions)

# Reusable HFold implementation. main.cc must not be included here because the
# unit-test executable gets its main() from GoogleTest.
add_library(HFoldLib STATIC
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

target_include_directories(HFoldLib PUBLIC
    "${PROJECT_SOURCE_DIR}/src"
)

target_compile_definitions(HFoldLib PUBLIC
    PARAMS_DIR="${PROJECT_SOURCE_DIR}/params"
)

target_link_libraries(HFoldLib PUBLIC
    RNA
    HFoldBuildOptions
)

# Command-line program.
add_executable(HFold src/main.cc)
target_link_libraries(HFold PRIVATE HFoldLib)

install(TARGETS HFold
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

if (NOT ENABLE_TESTS)
    message(STATUS "Unit tests are disabled. To enable them, set -DBUILD_TESTING=ON")
elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/googletest/CMakeLists.txt")
    enable_testing()

    set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)

    add_subdirectory(test/googletest EXCLUDE_FROM_ALL)
    add_subdirectory(test)
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


# Preserve the BUILD_PYTHON option observed in the original configuration.
# This recognizes the two conventional binding locations without affecting the
# normal C++ build when BUILD_PYTHON is OFF.
if(BUILD_PYTHON)
    if(EXISTS "${PROJECT_SOURCE_DIR}/python/CMakeLists.txt")
        add_subdirectory(python)
    elseif(EXISTS "${PROJECT_SOURCE_DIR}/bindings/python/CMakeLists.txt")
        add_subdirectory(bindings/python)
    else()
        message(FATAL_ERROR
            "BUILD_PYTHON is ON, but no Python binding CMakeLists.txt was found")
    endif()
endif()

# Source-formatting target. It always exists so `cmake --build build --target
# format` produces a useful result even when clang-format is not installed.
find_program(CLANG_FORMAT_EXECUTABLE NAMES clang-format)

file(GLOB_RECURSE HFOLD_FORMAT_FILES CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/src/*.c"
    "${PROJECT_SOURCE_DIR}/src/*.cc"
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
    "${PROJECT_SOURCE_DIR}/src/*.h"
    "${PROJECT_SOURCE_DIR}/src/*.hh"
    "${PROJECT_SOURCE_DIR}/src/*.hpp"
    "${PROJECT_SOURCE_DIR}/test/tests/*.c"
    "${PROJECT_SOURCE_DIR}/test/tests/*.cc"
    "${PROJECT_SOURCE_DIR}/test/tests/*.cpp"
    "${PROJECT_SOURCE_DIR}/test/tests/*.h"
    "${PROJECT_SOURCE_DIR}/test/tests/*.hh"
    "${PROJECT_SOURCE_DIR}/test/tests/*.hpp"
)

if(CLANG_FORMAT_EXECUTABLE)
    add_custom_target(format
        COMMAND "${CLANG_FORMAT_EXECUTABLE}" -i ${HFOLD_FORMAT_FILES}
        COMMENT "Formatting HFold source files"
        VERBATIM
    )
else()
    add_custom_target(format
        COMMAND "${CMAKE_COMMAND}" -E echo
                "clang-format was not found; install it to format the sources"
        VERBATIM
    )
endif()
