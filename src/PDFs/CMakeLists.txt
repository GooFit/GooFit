# This is only needed for IDE header display
file(GLOB GOOFIT_PDF_HEADERS
    "${PROJECT_SOURCE_DIR}/include/goofit/PDFs/*.h"
    "${PROJECT_SOURCE_DIR}/include/goofit/PDFs/detail/*.h"
    "${PROJECT_SOURCE_DIR}/include/goofit/PDFs/utilities/*.h"
    )

set(GOOPDF_FILE_LISTING
    detail/Globals.cu
    detail/Globals.cpp
    GooPdf.cpp
    GooPdf.cu
    MetricPointer.cu
    MetricTaker.cu
    ParameterContainer.cu
    PdfBase.cu
    utilities/VariableBinTransform1DPdf.cu
    ${GOOFIT_PDF_HEADERS}
    )


if(IS_NOT_CUDA OR NEW_CUDA)
    macro(goofit_add_pdf_library NAME)
        goofit_add_library(${NAME} ${ARGN})
        set_target_properties(${NAME} PROPERTIES FOLDER "pdfs")
        target_link_libraries(${NAME} PUBLIC PDFCore)
    endmacro()

    goofit_add_library(PDFCore ${GOOPDF_FILE_LISTING})
    target_link_libraries(PDFCore PUBLIC goofit_base)
    set_target_properties(PDFCore PROPERTIES FOLDER "pdfs")

    add_subdirectory(basic)
    add_subdirectory(combine)
    add_subdirectory(physics)

    add_library(PDFs INTERFACE)
    target_link_libraries(PDFs INTERFACE PDFBasic PDFCombine PDFPhysics)

# Workaround for FindCUDA (remove when CUDA 3.8+ is required)
else()
    macro(goofit_add_pdf_library NAME)
        set(GOOPDF_TEMP "${ARGN}")
        foreach(F IN LISTS GOOPDF_TEMP)
            if(IS_ABSOLUTE "${F}")
                list(APPEND GOOPDF_TEMP_2 "${F}")
            else()
                list(APPEND GOOPDF_TEMP_2 "${CMAKE_CURRENT_SOURCE_DIR}/${F}")
            endif()
        endforeach()

        set(GOOPDF_FILE_LISTING ${GOOPDF_FILE_LISTING} ${GOOPDF_TEMP_2} PARENT_SCOPE)
        unset(GOOPDF_TEMP)
        unset(GOOPDF_TEMP_2)
    endmacro()

    add_subdirectory(basic)
    add_subdirectory(combine)
    add_subdirectory(physics)

    goofit_add_library(PDFs ${GOOPDF_FILE_LISTING})
    target_link_libraries(PDFs PUBLIC goofit_base)
    set_target_properties(PDFs PROPERTIES FOLDER "pdfs")
endif()
