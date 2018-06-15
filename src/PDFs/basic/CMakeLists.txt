# This is only needed for IDE header display
file(GLOB_RECURSE GOOFIT_BASIC_HEADERS
    "${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/*.h"
    )

goofit_add_library(PDFBasic
    ArgusPdf.cu
    BWPdf.cu
    BifurGaussPdf.cu
    BinTransformPdf.cu
    CorrGaussianPdf.cu
    CrystalBallPdf.cu
    ExpGausPdf.cu
    ExpPdf.cu
    GaussianPdf.cu
    InterHistPdf.cu
    JohnsonSUPdf.cu
    KinLimitBWPdf.cu
    LandauPdf.cu
    NovosibirskPdf.cu
    PolynomialPdf.cu
    ScaledGaussianPdf.cu
    SmoothHistogramPdf.cu
    StepPdf.cu
    TrigThresholdPdf.cu
    VoigtianPdf.cu
    ${GOOFIT_BASIC_HEADERS}
    )

set_target_properties(PDFBasic PROPERTIES FOLDER "pdfs")
target_link_libraries(PDFBasic PUBLIC PDFCore)
