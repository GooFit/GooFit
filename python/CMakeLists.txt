pybind11_add_module(pygoofit pygoofit.cpp Variable.cpp BinnedDataSet.cpp)
target_link_libraries(pygoofit PUBLIC goofit_lib)
