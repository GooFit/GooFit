pybind11_add_module(goofit
    goofit.cpp
    Variable.cpp
    DataSet.cpp BinnedDataSet.cpp UnbinnedDataSet.cpp
    PdfBase.cpp
    PDFs/GooPdf.cpp
    PDFs/ExpPdf.cpp
    FitManager.cpp)

if(GOOFIT_DEVICE STREQUAL CUDA)
    target_include_directories(goofit PUBLIC ${CUDA_INCLUDE_DIRS})
endif()

target_link_libraries(goofit PUBLIC goofit_lib)
set_target_properties(goofit PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
