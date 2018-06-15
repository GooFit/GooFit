# This is only needed for IDE header display
file(GLOB_RECURSE GOOFIT_COMBINE_HEADERS
    "${PROJECT_SOURCE_DIR}/include/goofit/PDFs/combine/*.h"
    )

goofit_add_library(PDFCombine
    AddPdf.cu
    CompositePdf.cu
    ConvolutionPdf.cu
    EventWeightedAddPdf.cu
    MappedPdf.cu
    ProdPdf.cu
    ${GOOFIT_COMBINE_HEADERS}
)

set_target_properties(PDFCombine PROPERTIES FOLDER "pdfs")

target_link_libraries(PDFCombine PUBLIC PDFCore)
