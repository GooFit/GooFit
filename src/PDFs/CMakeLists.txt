
include_directories(
    ${PROJECT_SOURCE_DIR}/src/PDFs
    )

goofit_add_library(PDFs AllPdfs.cu)
target_link_libraries(PDFs goofit_base)

export(TARGETS PDFs APPEND FILE ${PROJECT_BINARY_DIR}/GooFitConfig.cmake)
