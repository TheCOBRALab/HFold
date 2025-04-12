cmake_minimum_required(VERSION 3.15)
project(HFold
    LANGUAGES C CXX
    DESCRIPTION "RNA folding prediction tool"
)

# Set standards before targets are created
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default to Release build if no configuration is specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Organize sources using modern list structure
set(hfold_sources
    src/main.cpp
    src/W_final.cpp 
    src/pseudo_loop.cpp 
    src/HFold.cpp
    src/cmdline.cpp
    src/Result.cpp	
    src/s_energy_matrix.cpp	
    src/Hotspot.cpp
    src/sparse_tree.cpp
)

set(vienna_sources
    # Constraints
    src/ViennaRNA/constraints/constraints.c
    src/ViennaRNA/constraints/hard.c
    src/ViennaRNA/constraints/SHAPE.c
    src/ViennaRNA/constraints/soft.c
    
    # Datastructures
    src/ViennaRNA/datastructures/basic_datastructures.c
    src/ViennaRNA/datastructures/char_stream.c
    src/ViennaRNA/datastructures/lists.c
    
    # I/O
    src/ViennaRNA/io/file_formats.c
    src/ViennaRNA/io/io_utils.c
    
    # Landscape
    src/ViennaRNA/landscape/move.c
    
    # Loops
    src/ViennaRNA/loops/external_pf.c
    src/ViennaRNA/loops/external.c
    src/ViennaRNA/loops/hairpin.c
    src/ViennaRNA/loops/internal.c
    src/ViennaRNA/loops/multibranch.c
    
    # Parameters
    src/ViennaRNA/params/default.c
    src/ViennaRNA/params/params.c
    src/ViennaRNA/params/io.c
    
    # Utilities
    src/ViennaRNA/utils/cpu.c
    src/ViennaRNA/utils/higher_order_functions.c
    src/ViennaRNA/utils/string_utils.c
    src/ViennaRNA/utils/structure_utils.c
    src/ViennaRNA/utils/utils.c
    
    # Core ViennaRNA
    src/ViennaRNA/alphabet.c
    src/ViennaRNA/boltzmann_sampling.c
    src/ViennaRNA/centroid.c
    src/ViennaRNA/cofold.c
    src/ViennaRNA/commands.c
    src/ViennaRNA/dp_matrices.c
    src/ViennaRNA/equilibrium_probs.c
    src/ViennaRNA/eval.c
    src/ViennaRNA/fold_compound.c
    src/ViennaRNA/fold.c
    src/ViennaRNA/gquad.c
    src/ViennaRNA/grammar.c
    src/ViennaRNA/MEA.c
    src/ViennaRNA/mfe.c
    src/ViennaRNA/mm.c
    src/ViennaRNA/model.c
    src/ViennaRNA/part_func.c
    src/ViennaRNA/ribo.c
    src/ViennaRNA/sequence.c
    src/ViennaRNA/subopt.c
    src/ViennaRNA/unstructured_domains.c
)

# RNA static library (ViennaRNA)
add_library(RNA STATIC ${vienna_sources})
target_include_directories(RNA PUBLIC 
    src
)

target_compile_options(RNA PRIVATE
    $<$<CONFIG:Debug>:-g -Wall -Wextra>
    $<$<CONFIG:Release>:-O3 -flto>
)

# HFold executable
add_executable(HFold ${hfold_sources})
target_link_libraries(HFold PRIVATE RNA)
target_compile_options(HFold PRIVATE
    $<$<CONFIG:Debug>:-g -Wall -Wextra> # cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Debug
    $<$<CONFIG:Release>:-O3 -flto>      # cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release
)

# Installation (For Conda package)
install(TARGETS HFold RUNTIME DESTINATION bin)