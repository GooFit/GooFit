# This is only needed for IDE header display
file(GLOB_RECURSE GOOFIT_COMBINE_HEADERS "${PROJECT_SOURCE_DIR}/include/goofit/PDFs/combine/*.h")

goofit_add_pdf_library(
  PDFCombine
  AddPdf.cu
  CompositePdf.cu
  ConvolutionPdf.cu
  EventWeightedAddPdf.cu
  MappedPdf.cu
  ProdPdf.cu
  ${GOOFIT_COMBINE_HEADERS})
