pybind11_add_module(goofit goofit.cpp
    Variable.cpp
    DataSet.cpp BinnedDataSet.cpp UnbinnedDataSet.cpp
    PdfBase.cpp
    PDFs/GooPdf.cpp
    PDFs/ExpPdf.cpp
    FitManager.cpp)
target_link_libraries(goofit PUBLIC goofit_lib)
set_target_properties(goofit PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
