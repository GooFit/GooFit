
set(DATA_FILE "${CMAKE_CURRENT_SOURCE_DIR}/zachFit_datafiles.tgz")
 
if(NOT EXISTS ${DATA_FILE})
    message(STATUS "Downloading zachFit data files")
    file(DOWNLOAD
        https://github.com/GooFit/GooFit/releases/download/v1.0.0/zachFit_dataFiles.tgz
        ${DATA_FILE}
        SHOW_PROGRESS
        EXPECTED_HASH SHA256=66340df3799392272ee32d6508323d12f8fcca7eb09264ef4f140ae5246f644e
        )
endif()

if(NOT EXISTS 
        "${CMAKE_CURRENT_SOURCE_DIR}/datafiles/dstwidth_kpi_resMC.dat"
        )
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E tar xzf ${DATA_FILE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
endif()


goofit_add_directory()
goofit_add_executable(zachFit zachFit.cu)
