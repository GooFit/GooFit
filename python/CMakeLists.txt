pybind11_add_module(goofit
    goofit.cpp
    goofit/Variable.cpp
    goofit/DataSet.cpp
    goofit/BinnedDataSet.cpp
    goofit/UnbinnedDataSet.cpp
    goofit/PdfBase.cpp
    PDFs/GooPdf.cpp
    PDFs/basic/ExpPdf.cpp
    goofit/FitManager.cpp)

if(GOOFIT_DEVICE STREQUAL CUDA)
    target_include_directories(goofit PUBLIC ${CUDA_INCLUDE_DIRS})
endif()

target_link_libraries(goofit PUBLIC goofit_lib)
set_target_properties(goofit PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set_target_properties(goofit PROPERTIES FOLDER python)

add_custom_target(pygoofit_tests ALL
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/tests ${PROJECT_BINARY_DIR}/pytests)
set_target_properties(pygoofit_tests PROPERTIES FOLDER python)
add_custom_target(pygoofit_examples ALL
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/examples ${PROJECT_BINARY_DIR}/pyexamples)
set_target_properties(pygoofit_examples PROPERTIES FOLDER python)
