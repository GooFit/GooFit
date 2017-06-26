pybind11_add_module(goofit
    goofit.cpp
    goofit/Variable.cpp
    goofit/DataSet.cpp
    goofit/BinnedDataSet.cpp
    goofit/UnbinnedDataSet.cpp
    goofit/PdfBase.cpp
    PDFs/GooPdf.cpp
    PDFs/basic/ArgusPdf.cpp
    PDFs/basic/BifurGaussPdf.cpp
    PDFs/basic/BinTransformPdf.cpp
    PDFs/basic/BWPdf.cpp
    PDFs/basic/CorrGaussianPdf.cpp
    PDFs/basic/CrystalBallPdf.cpp
    PDFs/basic/ExpGausPdf.cpp
    PDFs/basic/ExpPdf.cpp
    PDFs/basic/GaussianPdf.cpp
    PDFs/basic/InterHistPdf.cpp
    PDFs/basic/JohnsonSUPdf.cpp
    PDFs/basic/KinLimitBWPdf.cpp
    PDFs/basic/LandauPdf.cpp
    PDFs/basic/NovosibirskPdf.cpp
    PDFs/basic/PolynomialPdf.cpp
    PDFs/basic/ScaledGaussianPdf.cpp
    PDFs/basic/SmoothHistogramPdf.cpp
    PDFs/basic/StepPdf.cpp
    PDFs/basic/TrigThresholdPdf.cpp
    PDFs/basic/VoigtianPdf.cpp
    PDFs/combine/AddPdf.cpp
    PDFs/combine/CompositePdf.cpp
    PDFs/combine/ConvolutionPdf.cpp
    PDFs/combine/EventWeightedAddPdf.cpp
    PDFs/combine/MappedPdf.cpp
    PDFs/combine/ProdPdf.cpp
    PDFs/physics/DalitzPlotPdf.cpp
    PDFs/physics/DalitzVetoPdf.cpp
    PDFs/physics/DP4Pdf.cpp
    PDFs/physics/IncoherentSumPdf.cpp
    PDFs/physics/LineshapesPdf.cpp
    PDFs/physics/ResonancePdf.cpp
    PDFs/physics/Tddp4Pdf.cpp
    PDFs/physics/TddpPdf.cpp
    goofit/FitManager.cpp)

if(GOOFIT_DEVICE STREQUAL CUDA)
    target_include_directories(goofit PUBLIC ${CUDA_INCLUDE_DIRS})
endif()

target_link_libraries(goofit PUBLIC goofit_lib)
set_target_properties(goofit PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set_target_properties(goofit PROPERTIES FOLDER python)

add_custom_target(pygoofit_tests ALL
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/tests ${PROJECT_BINARY_DIR}/pytests)
set_target_properties(pygoofit_tests PROPERTIES FOLDER python)
add_custom_target(pygoofit_examples ALL
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/examples ${PROJECT_BINARY_DIR}/pyexamples)
set_target_properties(pygoofit_examples PROPERTIES FOLDER python)
