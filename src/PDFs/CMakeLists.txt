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
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/DalitzPlotter.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/SpinFactors.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/MixingTimeResolution_Aux.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/TruthResolution_Aux.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/ThreeGaussResolution_Aux.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/EvalVar.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/UserUtils.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/physics/SpinHelper.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/utility/VariableBinTransform1DPdf.h
    ${PROJECT_SOURCE_DIR}/include/goofit/PDFs/detail/ComplexUtils.h
)

if(GOOFIT_KMATRIX AND GOOFIT_DEVICE STREQUAL CUDA)
    list(APPEND GOOPDF_HEADERS ${PROJECT_SOURCE_DIR}/include/goofit/detail/compute_inverse5.h)
endif()

get_property(goo_extra GLOBAL PROPERTY GOOFIT_EXTRA_PDFS)

if(GOOFIT_SEPARATE_COMP)
    goofit_add_library(PDFs
        GooPdf.cu
        MetricTaker.cu
        PdfBase.cu
        basic/CrystalBallPdf.cu
        basic/BifurGaussPdf.cu
        physics/LineshapesPdf.cu
        basic/CorrGaussianPdf.cu
        basic/LandauPdf.cu
        physics/DP4Pdf.cu
        basic/GaussianPdf.cu
        combine/AddPdf.cu
        basic/VoigtianPdf.cu
        combine/MappedPdf.cu
        physics/Tddp4Pdf.cu
        physics/DalitzVetoPdf.cu
        basic/JohnsonSUPdf.cu
        combine/ProdPdf.cu
        basic/PolynomialPdf.cu
        physics/IncoherentSumPdf.cu
        combine/EventWeightedAddPdf.cu
        basic/ScaledGaussianPdf.cu
        basic/StepPdf.cu
        basic/KinLimitBWPdf.cu
        basic/NovosibirskPdf.cu
        physics/ResonancePdf.cu
        basic/BWPdf.cu
        basic/BinTransformPdf.cu
        basic/TrigThresholdPdf.cu
        basic/SmoothHistogramPdf.cu
        basic/ExpPdf.cu
        basic/InterHistPdf.cu
        physics/DalitzPlotPdf.cu
        physics/TddpPdf.cu
        basic/ArgusPdf.cu
        basic/ExpGausPdf.cu
        combine/CompositePdf.cu
        combine/ConvolutionPdf.cu
        physics/DalitzPlotHelpers.cu
        physics/SpinFactors.cu
        physics/MixingTimeResolution_Aux.cu
        physics/TruthResolution_Aux.cu
        physics/ThreeGaussResolution_Aux.cu
        utility/VariableBinTransform1DPdf.cu
        ${GOOPDF_HEADERS}
        ${goo_extra}
    )

else()
    if(goo_extra)
        message(FATAL_ERROR "You must use separable compilation (GOOFIT_SEPARATE_COMP=ON) to add PDFs!)")
    endif()

    target_include_directories(PDFs PUBLIC ${PROJECT_SOURCE_DIR}/src/PDFs)

    goofit_add_library(PDFs AllPdfs.cu ${GOOPDF_HEADERS})

endif()


target_link_libraries(PDFs PUBLIC goofit_base)
