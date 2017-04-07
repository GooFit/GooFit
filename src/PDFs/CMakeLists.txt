set(GOOPDF_HEADERS
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/GooPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/CrystalBallPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/BifurGaussPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/LineshapesPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/CorrGaussianPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/LandauPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/DP4Pdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/GaussianPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/AddPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/VoigtianPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/MappedPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/Tddp4Pdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/DalitzVetoPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/JohnsonSUPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/ProdPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/PolynomialPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/IncoherentSumPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/EventWeightedAddPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/ScaledGaussianPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/StepPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/KinLimitBWPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/NovosibirskPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/ResonancePdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/BWPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/BinTransformPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/TrigThresholdPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/SmoothHistogramPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/ExpPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/InterHistPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/DalitzPlotPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/TddpPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/ArgusPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/ExpGausPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/CompositePdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/ConvolutionPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/DalitzPlotHelpers.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/SpinFactors.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/MixingTimeResolution_Aux.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/TruthResolution_Aux.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/ThreeGaussResolution_Aux.h
)

get_property(goo_extra GLOBAL PROPERTY GOOFIT_EXTRA_PDFS)

if(GOOFIT_SEPARATE_COMP)
    goofit_add_library(PDFs
        GooPdf.cu
        CrystalBallPdf.cu
        BifurGaussPdf.cu
        LineshapesPdf.cu
        CorrGaussianPdf.cu
        LandauPdf.cu
        DP4Pdf.cu
        GaussianPdf.cu
        AddPdf.cu
        VoigtianPdf.cu
        MappedPdf.cu
        Tddp4Pdf.cu
        DalitzVetoPdf.cu
        JohnsonSUPdf.cu
        ProdPdf.cu
        PolynomialPdf.cu
        IncoherentSumPdf.cu
        EventWeightedAddPdf.cu
        ScaledGaussianPdf.cu
        StepPdf.cu
        KinLimitBWPdf.cu
        NovosibirskPdf.cu
        ResonancePdf.cu
        BWPdf.cu
        BinTransformPdf.cu
        TrigThresholdPdf.cu
        SmoothHistogramPdf.cu
        ExpPdf.cu
        InterHistPdf.cu
        DalitzPlotPdf.cu
        TddpPdf.cu
        ArgusPdf.cu
        ExpGausPdf.cu
        CompositePdf.cu
        ConvolutionPdf.cu
        DalitzPlotHelpers.cu
        SpinFactors.cu
        MixingTimeResolution_Aux.cu
        TruthResolution_Aux.cu
        ThreeGaussResolution_Aux.cu
        ${GOOPDF_HEADERS}
        ${goo_extra})
else()
    if(goo_extra)
        message(FATAL_ERROR "You must use separable compilation (GOOFIT_SEPARATE_COMP=ON) to add PDFs!)")
    endif()

    include_directories(${PROJECT_SOURCE_DIR}/src/PDFs)
    
    goofit_add_library(PDFs AllPdfs.cu ${GOOPDF_HEADERS})

endif()

target_link_libraries(PDFs goofit_base)
install(TARGETS PDFs EXPORT GooFitTargets DESTINATION lib)
