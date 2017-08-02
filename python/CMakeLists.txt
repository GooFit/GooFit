goofit_add_library(_goofit MODULE
    goofit.cpp
    goofit/Variable.cpp
    goofit/DataSet.cpp
    goofit/BinnedDataSet.cpp
    goofit/UnbinnedDataSet.cpp
    goofit/PdfBase.cpp
    goofit/Version.cpp
    PDFs/GooPdf.cu
    PDFs/basic/ArgusPdf.cu
    PDFs/basic/BifurGaussPdf.cu
    PDFs/basic/BinTransformPdf.cu
    PDFs/basic/BWPdf.cu
    PDFs/basic/CorrGaussianPdf.cu
    PDFs/basic/CrystalBallPdf.cu
    PDFs/basic/ExpGausPdf.cu
    PDFs/basic/ExpPdf.cu
    PDFs/basic/GaussianPdf.cu
    PDFs/basic/InterHistPdf.cu
    PDFs/basic/JohnsonSUPdf.cu
    PDFs/basic/KinLimitBWPdf.cu
    PDFs/basic/LandauPdf.cu
    PDFs/basic/NovosibirskPdf.cu
    PDFs/basic/PolynomialPdf.cu
    PDFs/basic/ScaledGaussianPdf.cu
    PDFs/basic/SmoothHistogramPdf.cu
    PDFs/basic/StepPdf.cu
    PDFs/basic/TrigThresholdPdf.cu
    PDFs/basic/VoigtianPdf.cu
    PDFs/combine/AddPdf.cu
    PDFs/combine/CompositePdf.cu
    PDFs/combine/ConvolutionPdf.cu
    PDFs/combine/EventWeightedAddPdf.cu
    PDFs/combine/MappedPdf.cu
    PDFs/combine/ProdPdf.cu
    PDFs/physics/DalitzPlotPdf.cu
    PDFs/physics/DalitzVetoPdf.cu
    PDFs/physics/DP4Pdf.cu
    PDFs/physics/IncoherentSumPdf.cu
    PDFs/physics/LineshapesPdf.cu
    PDFs/physics/MixingTimeResolution_Aux.cu
    PDFs/physics/ResonancePdf.cu
    PDFs/physics/SpinFactors.cu
    PDFs/physics/Tddp4Pdf.cu
    PDFs/physics/TddpPdf.cu
    PDFs/physics/ThreeGaussResolution_Aux.cu
    PDFs/physics/TruthResolution_Aux.cu
    goofit/FitManager.cpp)

target_link_libraries(_goofit pybind11::module)

set_target_properties(_goofit PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                         SUFFIX "${PYTHON_MODULE_EXTENSION}")

# Only need to add CUDA incdir if build with CUDA
# (Issue with pre-CMake 3.8 design)
if(GOOFIT_DEVICE STREQUAL CUDA)
    target_include_directories(_goofit PUBLIC ${CUDA_INCLUDE_DIRS})
endif()

target_link_libraries(_goofit goofit_lib)
set_target_properties(_goofit PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/goofit)
set_target_properties(_goofit PROPERTIES FOLDER python)

install(TARGETS _goofit LIBRARY DESTINATION goofit)


file(WRITE "${CMAKE_BINARY_DIR}/goofit/__init__.py" "from ._goofit import *")
install(FILES "${CMAKE_BINARY_DIR}/goofit/__init__.py" DESTINATION "goofit")

add_custom_target(pygoofit_tests ALL
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/tests ${PROJECT_BINARY_DIR}/pytests)
set_target_properties(pygoofit_tests PROPERTIES FOLDER python)
add_custom_target(pygoofit_examples ALL
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/examples ${PROJECT_BINARY_DIR}/pyexamples)
set_target_properties(pygoofit_examples PROPERTIES FOLDER python)
