# This is only needed for IDE header display
file(GLOB GOOFIT_PDF_HEADERS "${PROJECT_SOURCE_DIR}/include/goofit/PDFs/*.h"
     "${PROJECT_SOURCE_DIR}/include/goofit/PDFs/detail/*.h"
     "${PROJECT_SOURCE_DIR}/include/goofit/PDFs/utilities/*.h")

set(GOOPDF_FILE_LISTING
    detail/Globals.cu
    detail/Globals.cpp
    GooPdf.cpp
    GooPdf.cu
    MetricPointer.cu
    MetricTaker.cu
    ParameterContainer.cu
    PdfBase.cu
    utilities/VariableBinTransform1DPdf.cu
    utilities/DebugTools.cu
    ${GOOFIT_PDF_HEADERS})

macro(goofit_add_pdf_library NAME)
  goofit_add_library(${NAME} ${ARGN})
  set_target_properties(${NAME} PROPERTIES FOLDER "pdfs")
  target_link_libraries(${NAME} PUBLIC PDFCore)
endmacro()

goofit_add_library(PDFCore ${GOOPDF_FILE_LISTING})
target_link_libraries(PDFCore PUBLIC goofit_base)
set_target_properties(PDFCore PROPERTIES FOLDER "pdfs")

add_subdirectory(basic)
add_subdirectory(combine)
add_subdirectory(physics)

add_library(PDFs INTERFACE)
target_link_libraries(PDFs INTERFACE PDFBasic PDFCombine PDFPhysics)
