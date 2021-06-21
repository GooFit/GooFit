add_library(
  minuit2 MODULE
  main.cpp
  FCNBase.cpp
  FunctionMinimum.cpp
  MinuitParameter.cpp
  MnApplication.cpp
  MnMigrad.cpp
  MnPrint.cpp
  MnUserCovariance.cpp
  MnUserParameterState.cpp
  MnUserParameters.cpp
  MnScan.cpp
  MnMinos.cpp)

process_pybind_module(minuit2)
target_link_libraries(minuit2 PUBLIC ROOT::Minuit2)
target_link_libraries(minuit2 PRIVATE pybind11::module)
target_include_directories(minuit2 PRIVATE "${GOOFIT_BINARY_DIR}/include")
