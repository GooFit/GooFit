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

function(GOOFIT_ADD_PDF_DOCS PDFFILE)
    set(output_file "${GOOFIT_BINARY_DIR}/include/goofit/docs/${PDFFILE}")
    set(input_file "${GOOFIT_SOURCE_DIR}/include/goofit/${PDFFILE}")

    add_custom_command(
      OUTPUT
        "${output_file}"
      COMMAND
        ${CMAKE_COMMAND} -DOUTFILE=${output_file} -DINFILE=${input_file} -P "${CMAKE_CURRENT_LIST_DIR}/make_docs_header.cmake"
      VERBATIM
      MAIN_DEPENDENCY
        "${input_file}")

    target_sources(_goofit PRIVATE "${output_file}")
endfunction()

# Adding the Landau helper

add_library(landau MODULE extras/landau.cpp)
process_pybind_module(landau)
target_include_directories(landau PUBLIC ../include)

add_subdirectory(Minuit2)

goofit_add_library(_goofit MODULE
    include/goofit/PyProps.h
    include/goofit/Python.h
    goofit.cpp
    goofit/HelpPrinter.cpp
    goofit/Variable.cpp
    goofit/DataSet.cpp
    goofit/BinnedDataSet.cpp
    goofit/UnbinnedDataSet.cpp
    goofit/PdfBase.cpp
    goofit/Version.cpp
    goofit/FitControl.cpp
    goofit/Application.cpp
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
    PDFs/physics/DalitzPlotHelpers.cpp
    PDFs/physics/DalitzPlotter.cpp
    PDFs/physics/DP4Pdf.cu
    PDFs/physics/IncoherentSumPdf.cpp
    PDFs/physics/LineshapesPdf.cpp
    PDFs/physics/MixingTimeResolution_Aux.cpp
    PDFs/physics/ResonancePdf.cpp
    PDFs/physics/SpinFactors.cpp
    PDFs/physics/Tddp4Pdf.cu
    PDFs/physics/TddpPdf.cpp
    PDFs/physics/ThreeGaussResolution_Aux.cpp
    PDFs/physics/TruthResolution_Aux.cpp
    PDFs/utility/VariableBinTransform1DPdf.cpp
    goofit/FitManager.cpp)


goofit_add_pdf_docs(PDFs/basic/ArgusPdf.h)
goofit_add_pdf_docs(PDFs/basic/BifurGaussPdf.h)
goofit_add_pdf_docs(PDFs/basic/BWPdf.h)
goofit_add_pdf_docs(PDFs/basic/CorrGaussianPdf.h)
goofit_add_pdf_docs(PDFs/basic/CrystalBallPdf.h)
goofit_add_pdf_docs(PDFs/basic/ExpGausPdf.h)
goofit_add_pdf_docs(PDFs/basic/ExpPdf.h)
goofit_add_pdf_docs(PDFs/basic/GaussianPdf.h)
goofit_add_pdf_docs(PDFs/basic/InterHistPdf.h)
goofit_add_pdf_docs(PDFs/basic/JohnsonSUPdf.h)
goofit_add_pdf_docs(PDFs/basic/KinLimitBWPdf.h)
goofit_add_pdf_docs(PDFs/basic/LandauPdf.h)
goofit_add_pdf_docs(PDFs/basic/NovosibirskPdf.h)
goofit_add_pdf_docs(PDFs/basic/PolynomialPdf.h)
goofit_add_pdf_docs(PDFs/basic/ScaledGaussianPdf.h)
goofit_add_pdf_docs(PDFs/basic/SmoothHistogramPdf.h)
goofit_add_pdf_docs(PDFs/basic/StepPdf.h)
goofit_add_pdf_docs(PDFs/basic/VoigtianPdf.h)
goofit_add_pdf_docs(PDFs/basic/TrigThresholdPdf.h)
goofit_add_pdf_docs(PDFs/basic/BinTransformPdf.h)

goofit_add_pdf_docs(PDFs/combine/AddPdf.h)
goofit_add_pdf_docs(PDFs/combine/CompositePdf.h)
goofit_add_pdf_docs(PDFs/combine/ConvolutionPdf.h)
goofit_add_pdf_docs(PDFs/combine/EventWeightedAddPdf.h)
goofit_add_pdf_docs(PDFs/combine/MappedPdf.h)
goofit_add_pdf_docs(PDFs/combine/ProdPdf.h)

process_pybind_module(_goofit)
target_link_libraries(_goofit PUBLIC GooFit::GooFit)
target_include_directories(_goofit PUBLIC include)

# Only need to add CUDA incdir if built with FindCUDA

if(GOOFIT_DEVICE STREQUAL CUDA AND NOT NEW_CUDA)
    target_include_directories(_goofit PUBLIC ${CUDA_INCLUDE_DIRS})
endif()

# Create files in the python directory

file(WRITE "${CMAKE_BINARY_DIR}/goofit/__init__.py" [=[
from ._goofit import *
try:
    from IPython import get_ipython
    ip = get_ipython()
except ImportError:
    ip = None

if ip is not None:
    html_formatter = ip.display_formatter.formatters['text/markdown']
    html_formatter.for_type(type(PdfBase), lambda x: repr(x.help()) if hasattr(x, 'help') else None)
]=])
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
