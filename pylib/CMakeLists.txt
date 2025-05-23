
cmake_minimum_required(VERSION 3.15)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(hfold_python
    ${CMAKE_SOURCE_DIR}/src/W_final.cpp
    ${CMAKE_SOURCE_DIR}/src/pseudo_loop.cpp
    ${CMAKE_SOURCE_DIR}/src/HFold.cpp
    ${CMAKE_SOURCE_DIR}/src/cmdline.cpp
    ${CMAKE_SOURCE_DIR}/src/Result.cpp
    ${CMAKE_SOURCE_DIR}/src/s_energy_matrix.cpp
    ${CMAKE_SOURCE_DIR}/src/Hotspot.cpp
    ${CMAKE_SOURCE_DIR}/src/sparse_tree.cpp
    ${CMAKE_SOURCE_DIR}/pylib/bindings/pybind_module.cpp
)

target_include_directories(hfold_python
    PRIVATE
        ${VIENNARNA_INCLUDE_DIRECTORY}
        ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(hfold_python
    PRIVATE
        HFoldLib
        ${VIENNARNA_LIBRARY}
)

set_target_properties(hfold_python PROPERTIES
    OUTPUT_NAME hfold
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/pylib"
)

include(GNUInstallDirs)
install(TARGETS hfold_python
        LIBRARY DESTINATION hfold)
