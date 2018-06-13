add_subdirectory(basic)

add_subdirectory(combine)

add_subdirectory(physics)

goofit_add_library(PDFs
    detail/Globals.cu
    detail/Globals.cpp
    MetricPointer.cu
    GooPdf.cu
    GooPdf.cpp
    MetricTaker.cu
    ParameterContainer.cu
    PdfBase.cu
    utility/VariableBinTransform1DPdf.cu
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/detail/Globals.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/GooPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/CombinePdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/SimplePdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/ComplexComponent.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/MetricPointer.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/MetricTaker.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/ParameterContainer.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/detail/ComplexUtils.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/utility/VariableBinTransform1DPdf.h
)

target_link_libraries(PDFs PUBLIC PdfBasic PdfCombine PdfPhysics goofit_base)
