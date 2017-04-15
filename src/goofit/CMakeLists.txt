
goofit_add_library(Variable
    Variable.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/Variable.h
    )

goofit_add_library(FitManager
    FitManagerMinuit1.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/FitManager.h
    ${PROJECT_SOURCE_DIR}/include/goofit/FitManagerMinuit1.h
    )

goofit_add_library(FCN
    FCN.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/FCN.h
    )

goofit_add_library(FitManager2
    FitManagerMinuit2.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/FitManagerMinuit2.h
    )

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
    FitManager FunctorWriter FitControl
    PdfBase FitManager2 FCN
    )

add_library(goofit_base INTERFACE)
target_link_libraries(goofit_base INTERFACE ${lib_list})
install(TARGETS goofit_base ${lib_list} EXPORT GooFitTargets DESTINATION lib)
