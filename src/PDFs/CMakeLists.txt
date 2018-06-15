# This is only needed for IDE header display
file(GLOB GOOFIT_PDF_HEADERS
    "${PROJECT_SOURCE_DIR}/include/goofit/PDFs/*.h"
    "${PROJECT_SOURCE_DIR}/include/goofit/PDFs/detail/*.h"
    "${PROJECT_SOURCE_DIR}/include/goofit/PDFs/utilities/*.h"
    )

goofit_add_library(PDFCore
    detail/Globals.cu
    detail/Globals.cpp
    GooPdf.cpp
    GooPdf.cu
    MetricPointer.cu
    MetricTaker.cu
    ParameterContainer.cu
    PdfBase.cu
    utilities/VariableBinTransform1DPdf.cu
    ${GOOFIT_PDF_HEADERS}
)
target_link_libraries(PDFCore PUBLIC goofit_base)
set_target_properties(PDFCore PROPERTIES FOLDER "pdfs")

add_subdirectory(basic)
add_subdirectory(combine)
add_subdirectory(physics)

add_library(PDFs INTERFACE)
target_link_libraries(PDFs INTERFACE PDFBasic PDFCombine PDFPhysics)
