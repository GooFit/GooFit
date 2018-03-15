# Function to make adding a pybind_module easy and language agnostic
function(process_pybind_module LNAME)
    target_link_libraries(${LNAME} PUBLIC pybind11::module)

    set_target_properties(${LNAME} PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"
                                             SUFFIX "${PYTHON_MODULE_EXTENSION}")

    set_target_properties(${LNAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/goofit)
    set_target_properties(${LNAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY_DEBUG ${PROJECT_BINARY_DIR}/goofit)
    set_target_properties(${LNAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY_RELEASE ${PROJECT_BINARY_DIR}/goofit)
    set_target_properties(${LNAME} PROPERTIES FOLDER python)

    install(TARGETS ${LNAME} LIBRARY DESTINATION goofit)
endfunction()

# Adding the Landau helper

add_library(landau MODULE extras/landau.cpp)
process_pybind_module(landau)
target_include_directories(landau PUBLIC ../include)

add_subdirectory(Minuit2)

goofit_add_library(_goofit MODULE
    goofit.cpp
    Minuit2/FunctionMinimum.cpp
    goofit/props.h
    goofit/Variable.cpp
    goofit/DataSet.cpp
    goofit/BinnedDataSet.cpp
    goofit/UnbinnedDataSet.cpp
    goofit/PdfBase.cpp
    goofit/Version.cpp
    goofit/FitControl.cpp
    goofit/Application.cpp
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
    PDFs/physics/DalitzPlotHelpers.cu
    PDFs/physics/DalitzPlotter.cu
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
    PDFs/utility/VariableBinTransform1DPdf.cu
    goofit/FitManager.cpp)

process_pybind_module(_goofit)
target_link_libraries(_goofit PUBLIC GooFit::GooFit)

# Only need to add CUDA incdir if built with FindCUDA

if(GOOFIT_DEVICE STREQUAL CUDA AND NOT NEW_CUDA)
    target_include_directories(_goofit PUBLIC ${CUDA_INCLUDE_DIRS})
endif()

# Create files in the python directory

file(WRITE "${CMAKE_BINARY_DIR}/goofit/__init__.py" "from ._goofit import *")
install(FILES "${CMAKE_BINARY_DIR}/goofit/__init__.py" DESTINATION "goofit")

file(WRITE "${CMAKE_BINARY_DIR}/goofit/__main__.py" "from ._goofit import print_goofit_info; print_goofit_info()")
install(FILES "${CMAKE_BINARY_DIR}/goofit/__main__.py" DESTINATION "goofit")

# Adding links to files for simplicity of running

add_custom_target(pygoofit_tests ALL
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/tests ${PROJECT_BINARY_DIR}/pytests)
set_target_properties(pygoofit_tests PROPERTIES FOLDER python)
add_custom_target(pygoofit_examples ALL
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/examples ${PROJECT_BINARY_DIR}/pyexamples)
set_target_properties(pygoofit_examples PROPERTIES FOLDER python)

if(NOT SKBUILD)
    # Add a message about the PYTHON_PATH
    add_custom_command(TARGET _goofit POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan
                       "Python bindings built, use: export PYTHONPATH=${PROJECT_BINARY_DIR}")
endif()

