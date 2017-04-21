
goofit_add_library(Variable
    Variable.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/Variable.h
    )

if(ROOT_FOUND)
goofit_add_library(FitManager1
    FitManagerMinuit1.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/fitting/FitManagerMinuit1.h
    )
endif()

goofit_add_library(FCN
    Params.cc
    FCN.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/fitting/FCN.h
    ${PROJECT_SOURCE_DIR}/include/goofit/fitting/Params.h
    )

goofit_add_library(FitManager2
    FitManagerMinuit2.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/FitManager.h
    ${PROJECT_SOURCE_DIR}/include/goofit/fitting/FitManagerMinuit2.h
    )

target_link_libraries(FitManager2 rang CLI11)

goofit_add_library(Faddeeva
    Faddeeva.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/Faddeeva.h
    )

goofit_add_library(DataSet
    DataSet.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/DataSet.h
    )

goofit_add_library(BinnedDataSet
    BinnedDataSet.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/BinnedDataSet.h
    )

goofit_add_library(UnbinnedDataSet
    UnbinnedDataSet.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/UnbinnedDataSet.h
    )

goofit_add_library(FunctorWriter
    FunctorWriter.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/FunctorWriter.h
    )

goofit_add_library(FitControl
    FitControl.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/FitControl.h
    )


goofit_add_library(PdfBase
    PdfBase.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/PdfBase.h
    )


set(lib_list
    Variable Faddeeva DataSet
    BinnedDataSet UnbinnedDataSet
    FunctorWriter FitControl
    PdfBase FitManager2 FCN
    )

if(ROOT_FOUND)
    list(APPEND lib_list FitManager1)
endif()

add_library(goofit_base INTERFACE)
target_link_libraries(goofit_base INTERFACE ${lib_list})
