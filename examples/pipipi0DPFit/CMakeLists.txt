

set(DATA_FILE "${CMAKE_CURRENT_SOURCE_DIR}/alltxtfiles4pipipi0DPFit.tgz")

if(NOT EXISTS ${DATA_FILE})
file(DOWNLOAD
    https://github.com/GooFit/GooFit/releases/download/v1.0.0/alltxtfiles4pipipi0DPFit.tgz
    ${DATA_FILE}
    SHOW_PROGRESS
    EXPECTED_HASH SHA256=7c1a2ad656772d31fee9f6a43590e47574c9bdc32afce820820f30efca13803d
    )
endif()

if(NOT EXISTS 
        "${CMAKE_CURRENT_SOURCE_DIR}/dataFiles/toyPipipi0/dalitz_toyMC_000.txt"
        )
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E tar xzf ${DATA_FILE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
endif()


goofit_add_executable(pipipi0DPFit pipipi0DPFit.cu)

goofit_add_link(
alltxtfiles4pipipi0DPFit.tgz
bkg_2_pdf_sideband_6slices.txt
bkg_3_pdf.txt
bkg_4_pdf.txt
CMakeLists.txt
dataFiles
Makefile
pipipi0DPFit.cu
README
signal_sigma_2slices_pdf.txt
signal_sigma_3slices_pdf.txt
signal_sigma_4slices_pdf.txt
signal_sigma_5slices_pdf.txt
signal_sigma_6slices_pdf.txt
signal_sigma_7slices_pdf.txt
signal_sigma_8slices_pdf.txt)
