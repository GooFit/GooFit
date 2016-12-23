
set(goofit_pdfs_dirs
    ${PROJECT_SOURCE_DIR}/include/goofit
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs
    ${PROJECT_SOURCE_DIR}/MCBooster
    ${THRUST_INCLUDE_DIR}
    )

include_directories(${goofit_pdfs_dirs}
    ${PROJECT_SOURCE_DIR}/src
    )

goofit_add_library(PDFs AllPDFs.cu)

target_include_directories(PDFs PUBLIC ${GOOFIT_CUDA_OPTIONS})

export(TARGETS PDFs APPEND FILE ${PROJECT_BINARY_DIR}/GooFitConfig.cmake)
