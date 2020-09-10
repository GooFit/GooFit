set(DATA_FILE "${CMAKE_CURRENT_SOURCE_DIR}/alltxtfiles4pipipi0DPFit.tgz")

if(NOT EXISTS ${DATA_FILE})
  message(STATUS "Downloading pipipi0 data files")
  file(
    DOWNLOAD https://github.com/GooFit/GooFit/releases/download/v1.0.0/alltxtfiles4pipipi0DPFit.tgz
    ${DATA_FILE}
    EXPECTED_HASH SHA256=7c1a2ad656772d31fee9f6a43590e47574c9bdc32afce820820f30efca13803d
    SHOW_PROGRESS)
endif()

if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/signal_sigma_2slices_pdf.txt")
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${DATA_FILE}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

goofit_add_directory()
goofit_add_example(
  pipipi0DPFit
  TIME
  32
  ARGS
  toy
  -m
  10000)
