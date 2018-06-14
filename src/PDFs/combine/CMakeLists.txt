goofit_add_library(PdfCombine
    AddPdf.cu
    CompositePdf.cu
    ConvolutionPdf.cu
    EventWeightedAddPdf.cu
    MappedPdf.cu
    ProdPdf.cu
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/combine/AddPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/combine/CompositePdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/combine/ConvolutionPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/combine/EventWeightedAddPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/combine/MappedPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/combine/ProdPdf.h
)

target_link_libraries(PdfCombine PUBLIC PDFCore)
