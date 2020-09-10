add_library(
  _Core STATIC
  HelpPrinter.cpp
  Variable.cpp
  DataSet.cpp
  BinnedDataSet.cpp
  UnbinnedDataSet.cpp
  PdfBase.cpp
  Version.cpp
  FitControl.cpp
  Application.cpp
  FitManager.cpp)

target_link_libraries(_Core PRIVATE _goofit_python)
set_target_properties(_Core PROPERTIES FOLDER python)
