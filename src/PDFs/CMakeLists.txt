include_directories(
    ${PROJECT_SOURCE_DIR}/src/PDFs
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
        ${goo_extra})
else()
    if(goo_extra)
        message(FATAL_ERROR "You must use separable compilation (GOOFIT_SEPARATE_COMP=ON) to add PDFs!)")
    endif()
    goofit_add_library(PDFs AllPdfs.cu)
endif()

target_link_libraries(PDFs goofit_base)
install(TARGETS PDFs EXPORT GooFitTargets DESTINATION lib)
