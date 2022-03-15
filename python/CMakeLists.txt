# Function to make adding a pybind_module easy and language agnostic
function(process_pybind_module LNAME)
  set_target_properties(${LNAME} PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                            SUFFIX "${PYTHON_MODULE_EXTENSION}")

  set_target_properties(${LNAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/goofit)
  set_target_properties(${LNAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY_DEBUG
                                            ${PROJECT_BINARY_DIR}/goofit)
  set_target_properties(${LNAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY_RELEASE
                                            ${PROJECT_BINARY_DIR}/goofit)
  set_target_properties(${LNAME} PROPERTIES FOLDER python)

  install(TARGETS ${LNAME} LIBRARY DESTINATION goofit)
endfunction()

function(GOOFIT_ADD_PYTHON_DOCS TARGET PDFFILE)
  set(output_file "${GOOFIT_BINARY_DIR}/include/goofit/docs/${PDFFILE}")
  set(input_file "${GOOFIT_SOURCE_DIR}/include/goofit/${PDFFILE}")
  add_custom_command(
    OUTPUT "${output_file}"
    COMMAND ${CMAKE_COMMAND} -DOUTFILE=${output_file} -DINFILE=${input_file} -P
            "${GOOFIT_SOURCE_DIR}/python/make_docs_header.cmake"
    VERBATIM
    MAIN_DEPENDENCY "${input_file}")
  target_sources(${TARGET} PRIVATE "${output_file}")
endfunction()

# Adding the Landau helper

add_library(landau MODULE extras/landau.cpp)
process_pybind_module(landau)
target_link_libraries(landau PRIVATE pybind11::module)
target_include_directories(landau PUBLIC ../include)

# Adding Minuit2

add_subdirectory(Minuit2)

# Link to this library to get common includes

add_library(_goofit_python INTERFACE)
target_link_libraries(_goofit_python INTERFACE GooFit::GooFit pybind11::module)
target_include_directories(_goofit_python INTERFACE include)

# Adding GooFit Python bindings

goofit_add_library(_goofit MODULE include/goofit/PyProps.h include/goofit/Python.h goofit.cpp)

add_subdirectory(goofit)
add_subdirectory(PDFs)

process_pybind_module(_goofit)
target_link_libraries(_goofit PRIVATE _goofit_python _Core _Basic _Combine _Physics)

# Create files in the python directory

file(WRITE "${CMAKE_BINARY_DIR}/goofit/__init__.py" "from ._goofit import *")
install(FILES "${CMAKE_BINARY_DIR}/goofit/__init__.py" DESTINATION "goofit")

file(
  WRITE "${CMAKE_BINARY_DIR}/goofit/__main__.py"
  [=[
from __future__ import print_function
from ._goofit import goofit_info
print(goofit_info())
]=])

install(FILES "${CMAKE_BINARY_DIR}/goofit/__main__.py" DESTINATION "goofit")

# Adding links to files for simplicity of running

add_custom_target(
  pygoofit_tests ALL COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/tests
                             ${PROJECT_BINARY_DIR}/pytests)
set_target_properties(pygoofit_tests PROPERTIES FOLDER python)
add_custom_target(
  pygoofit_examples ALL
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/examples
          ${PROJECT_BINARY_DIR}/pyexamples)
set_target_properties(pygoofit_examples PROPERTIES FOLDER python)

# Add a message about the PYTHONPATH (if not built with pip)

if(NOT SKBUILD)
  add_custom_command(
    TARGET _goofit
    POST_BUILD
    COMMAND
      ${CMAKE_COMMAND} -E cmake_echo_color --cyan
      "Python bindings for ${PYTHON_EXECUTABLE} built, use: export PYTHONPATH=${PROJECT_BINARY_DIR}"
  )
endif()
