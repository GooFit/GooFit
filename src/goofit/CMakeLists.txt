
if(ROOT_FOUND)
goofit_add_library(FitManager1
    FitManagerMinuit1.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/fitting/FitManagerMinuit1.h
    )
endif()

goofit_add_library(FitManager2
    Params.cc
    FCN.cc
    FitManagerMinuit2.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/FitManager.h
    ${PROJECT_SOURCE_DIR}/include/goofit/fitting/FitManagerMinuit2.h
    ${PROJECT_SOURCE_DIR}/include/goofit/fitting/FCN.h
    ${PROJECT_SOURCE_DIR}/include/goofit/fitting/Params.h
    )

goofit_add_library(Faddeeva
    Faddeeva.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/Faddeeva.h
    )

goofit_add_library(DataSet
    DataSet.cc
    BinnedDataSet.cc
    UnbinnedDataSet.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/DataSet.h
    ${PROJECT_SOURCE_DIR}/include/goofit/BinnedDataSet.h
    ${PROJECT_SOURCE_DIR}/include/goofit/UnbinnedDataSet.h
    )


goofit_add_library(Variable
    Variable.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/Variable.h
    )

goofit_add_library(FunctorWriter
    FunctorWriter.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/FunctorWriter.h
    )


goofit_add_library(PdfBase
    PdfBase.cc
    Abort.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/detail/Abort.h
    ${PROJECT_SOURCE_DIR}/include/goofit/detail/ThrustOverride.h
    ${PROJECT_BINARY_DIR}/include/goofit/detail/ThrustOverrideConfig.h
    ${PROJECT_BINARY_DIR}/include/goofit/detail/Backtrace.h
    ${PROJECT_SOURCE_DIR}/include/goofit/Error.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PdfBase.h
    ${PROJECT_SOURCE_DIR}/include/goofit/FitControl.h
    ${PROJECT_SOURCE_DIR}/include/goofit/Application.h
    ${PROJECT_SOURCE_DIR}/include/goofit/Color.h
    ${PROJECT_SOURCE_DIR}/include/goofit/Log.h
    ${PROJECT_SOURCE_DIR}/include/goofit/GlobalCudaDefines.h
    )

target_link_libraries(PdfBase backtrace)


# Note: order is important on Linux, derived classes first
set(lib_list
    Faddeeva
    FunctorWriter
    PdfBase
    FitManager2
    DataSet
    Variable
    )

if(ROOT_FOUND)
    list(INSERT lib_list 0 FitManager1)
endif()

add_library(goofit_base INTERFACE)
target_link_libraries(goofit_base INTERFACE ${lib_list})
