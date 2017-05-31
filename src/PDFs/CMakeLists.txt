set(GOOPDF_HEADERS
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/GooPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/MetricTaker.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/CrystalBallPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/BifurGaussPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/LineshapesPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/CorrGaussianPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/LandauPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/DP4Pdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/GaussianPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/combine/AddPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/VoigtianPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/combine/MappedPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/Tddp4Pdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/DalitzVetoPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/JohnsonSUPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/combine/ProdPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/PolynomialPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/IncoherentSumPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/combine/EventWeightedAddPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/ScaledGaussianPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/StepPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/KinLimitBWPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/NovosibirskPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/ResonancePdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/BWPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/BinTransformPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/TrigThresholdPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/SmoothHistogramPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/ExpPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/InterHistPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/DalitzPlotPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/TddpPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/ArgusPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/basic/ExpGausPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/combine/CompositePdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/combine/ConvolutionPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/DalitzPlotHelpers.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/SpinFactors.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/MixingTimeResolution_Aux.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/TruthResolution_Aux.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/ThreeGaussResolution_Aux.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/EvalVar.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/SpinHelper.h
)

get_property(goo_extra GLOBAL PROPERTY GOOFIT_EXTRA_PDFS)

if(GOOFIT_SEPARATE_COMP)
    goofit_add_library(PDFs
        GooPdf.cu
        MetricTaker.cu
        PdfBase.cu
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
