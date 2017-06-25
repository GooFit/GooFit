pybind11_add_module(_goofit
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
    PDFs/physics/SpinFactors.cpp
    PDFs/physics/Tddp4Pdf.cpp
    PDFs/physics/TddpPdf.cpp
    goofit/FitManager.cpp)

# Only need to add CUDA incdir if build with CUDA
# (Issue with pre-CMake 3.8 design)
if(GOOFIT_DEVICE STREQUAL CUDA)
    target_include_directories(_goofit PUBLIC ${CUDA_INCLUDE_DIRS})
endif()

target_link_libraries(_goofit PUBLIC goofit_lib)
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
