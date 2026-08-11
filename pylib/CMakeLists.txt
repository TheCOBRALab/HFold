
cmake_minimum_required(VERSION 3.15)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(hfold_python
    ${CMAKE_SOURCE_DIR}/src/pseudo_loop.cc
    ${CMAKE_SOURCE_DIR}/src/pseudo_loop_trace.cc
    ${CMAKE_SOURCE_DIR}/src/HFold.cc
    ${CMAKE_SOURCE_DIR}/src/cmdline.cc
    ${CMAKE_SOURCE_DIR}/src/Result.cc
    ${CMAKE_SOURCE_DIR}/src/s_energy_matrix.cc
    ${CMAKE_SOURCE_DIR}/src/Hotspot.cc
    ${CMAKE_SOURCE_DIR}/src/sparse_tree.cc
    ${CMAKE_SOURCE_DIR}/src/SHAPE.cc
    ${CMAKE_SOURCE_DIR}/pylib/bindings/pybind_module.cc
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
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/pylib/hfold"
)

add_custom_command(TARGET hfold_python POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
    "$<TARGET_FILE:hfold_python>"
    "${CMAKE_SOURCE_DIR}/pylib/$<TARGET_FILE_NAME:hfold_python>"
)

include(GNUInstallDirs)
install(TARGETS hfold_python
        LIBRARY DESTINATION hfold)
