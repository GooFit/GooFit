
set(goofit_base_dirs 
    ${PROJECT_SOURCE_DIR}/include/goofit
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs
    ${PROJECT_SOURCE_DIR}/include/goofit/rootstuff)

include_directories(${goofit_base_dirs})

goofit_add_library(Variable
    Variable.cc
    )

goofit_add_library(FitManager
    FitManager.cc
    )

goofit_add_library(Faddeeva
    Faddeeva.cc
    )

goofit_add_library(DataSet
    DataSet.cc
    )

goofit_add_library(BinnedDataSet
    BinnedDataSet.cc
    )

goofit_add_library(UnbinnedDataSet
    UnbinnedDataSet.cc
    )

goofit_add_library(FunctorWriter
    FunctorWriter.cc
    )

goofit_add_library(FitControl
    FitControl.cc
    )


goofit_add_library(PdfBase
    PdfBase.cc
    )


set(lib_list
    Variable Faddeeva DataSet
    BinnedDataSet UnbinnedDataSet
    FitManager FunctorWriter FitControl
    PdfBase
    )

add_library(goofit_base STATIC)
target_link_libraries(goofit_base PUBLIC ${lib_list})
target_include_directories(goofit_base PUBLIC ${goofit_base_dirs})

export(TARGETS goofit_base ${lib_list} FILE ${PROJECT_BINARY_DIR}/GooFitConfig.cmake)
