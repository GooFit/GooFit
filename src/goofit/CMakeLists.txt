file(GLOB GOOFIT_EXTRA_HEADERS "${PROJECT_SOURCE_DIR}/include/goofit/detail/*.h"
     "${PROJECT_SOURCE_DIR}/include/goofit/utilities/*.h")

if(ROOT_FOUND)
  goofit_add_library(FitManager1 FitManagerMinuit1.cpp
                     ${PROJECT_SOURCE_DIR}/include/goofit/fitting/FitManagerMinuit1.h)
endif()

goofit_add_library(
  FitManager2
  Params.cpp
  FCN.cpp
  FitManagerMinuit2.cpp
  ${PROJECT_SOURCE_DIR}/include/goofit/FitManager.h
  ${PROJECT_SOURCE_DIR}/include/goofit/fitting/FitManagerMinuit2.h
  ${PROJECT_SOURCE_DIR}/include/goofit/fitting/FCN.h
  ${PROJECT_SOURCE_DIR}/include/goofit/fitting/Params.h)

goofit_add_library(Faddeeva Faddeeva.cpp ${PROJECT_SOURCE_DIR}/include/goofit/Faddeeva.h)

goofit_add_library(
  DataSet
  DataSet.cpp
  BinnedDataSet.cpp
  UnbinnedDataSet.cpp
  ${PROJECT_SOURCE_DIR}/include/goofit/DataSet.h
  ${PROJECT_SOURCE_DIR}/include/goofit/BinnedDataSet.h
  ${PROJECT_SOURCE_DIR}/include/goofit/UnbinnedDataSet.h)

goofit_add_library(Variable Variable.cpp ${PROJECT_SOURCE_DIR}/include/goofit/Variable.h)

goofit_add_library(FunctorWriter FunctorWriter.cpp
                   ${PROJECT_SOURCE_DIR}/include/goofit/FunctorWriter.h)

goofit_add_library(MathUtils MathUtils.cpp ${PROJECT_SOURCE_DIR}/include/goofit/MathUtils.h)

goofit_add_library(
  Application
  Application.cpp
  ${PROJECT_SOURCE_DIR}/include/goofit/Application.h
  ${PROJECT_SOURCE_DIR}/include/goofit/Color.h
  ${PROJECT_SOURCE_DIR}/include/goofit/Eigen.h
  ${PROJECT_SOURCE_DIR}/include/goofit/Error.h
  ${PROJECT_SOURCE_DIR}/include/goofit/FitControl.h
  ${PROJECT_SOURCE_DIR}/include/goofit/GlobalCudaDefines.h
  ${PROJECT_SOURCE_DIR}/include/goofit/Log.h
  ${GOOFIT_EXTRA_HEADERS})

if(GOOFIT_DEVICE STREQUAL CUDA)
  target_link_libraries(Application PRIVATE CUDA::Libs)
endif()

goofit_add_library(PdfBase PdfBase.cpp Abort.cpp ${PROJECT_SOURCE_DIR}/include/goofit/PdfBase.h
                   ${PROJECT_SOURCE_DIR}/include/goofit/detail/Abort.h)

target_link_libraries(PdfBase PUBLIC backtrace)
# Safe to use ROOT::Minuit2 even if ROOT not found
target_link_libraries(PdfBase PUBLIC ROOT::Minuit2)

# Note: order is important on Linux, derived classes first
set(lib_list
    Application
    Faddeeva
    FunctorWriter
    PdfBase
    FitManager2
    DataSet
    Variable
    MathUtils)

if(ROOT_FOUND)
  list(INSERT lib_list 0 FitManager1)
endif()

add_library(goofit_base INTERFACE)
target_link_libraries(goofit_base INTERFACE ${lib_list})
