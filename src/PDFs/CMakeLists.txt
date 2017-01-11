
include_directories(
    ${PROJECT_SOURCE_DIR}/src/PDFs
    )

goofit_add_library(PDFs AllPdfs.cu)
target_link_libraries(PDFs goofit_base)
install(TARGETS PDFs EXPORT GooFitTargets DESTINATION lib)
