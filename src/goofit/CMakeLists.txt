
goofit_add_library(Variable
    Variable.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/Variable.h
    )

goofit_add_library(FitManager
    FitManager.cc
    ${PROJECT_SOURCE_DIR}/include/goofit/FitManager.h
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
    PdfBase
    )

add_library(goofit_base INTERFACE)
target_link_libraries(goofit_base INTERFACE ${lib_list})
install(TARGETS goofit_base ${lib_list} EXPORT GooFitTargets DESTINATION lib)
