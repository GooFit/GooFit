# This macro sets up a goofit_add_example
# First argument is the name of the target, and it is followed by the source files
# If no source files present, assumes <name>.cu then <name>.cpp is a source file
#
# The options are:
#   NO_TEST:            Do not add a test
#   TEST_NAME <value>:  Add a test with this name instead of *_Examples
#   TIME <value>:       Add a time value to the test (integer, max 500)
#   ARGUMENTS <values>: Add arguments to the test call
#   EXISTING:           Don't add an example, just the test
#
function(GOOFIT_ADD_EXAMPLE EX_NAME)
  cmake_parse_arguments(GOOEX "EXISTING;NO_TEST" "TEST_NAME;TIME" "ARGS" ${ARGN})

  if(NOT GOOEX_EXISTING)
    if(GOOEX_UNPARSED_ARGUMENTS)
      goofit_add_executable(${EX_NAME} ${GOOEX_UNPARSED_ARGUMENTS})
    else()
      if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${EX_NAME}.cu")
        goofit_add_executable(${EX_NAME} ${EX_NAME}.cu)
      else()
        goofit_add_executable(${EX_NAME} ${EX_NAME}.cpp)
      endif()
    endif()
    set_target_properties(${EX_NAME} PROPERTIES FOLDER examples)
  endif()

  if(NOT GOOEX_NO_TEST)
    if(NOT GOOEX_TEST_NAME)
      set(GOOEX_TEST_NAME ${EX_NAME}_Example)
    endif()

    if(GOOEX_TEST_NAME)
      if(GOOEX_ARGS)
        if(GOOFIT_MPI)
          add_test(NAME ${GOOEX_TEST_NAME} COMMAND mpirun -np 2 ${EX_NAME} --nosplash
                                                   ${GOOEX_ARGS})
        else()
          add_test(NAME ${GOOEX_TEST_NAME} COMMAND ${EX_NAME} --nosplash ${GOOEX_ARGS})
        endif()
      else()
        if(GOOFIT_MPI)
          add_test(NAME ${GOOEX_TEST_NAME} COMMAND mpirun -np 2 ${EX_NAME} --nosplash)
        else()
          add_test(NAME ${GOOEX_TEST_NAME} COMMAND ${EX_NAME} --nosplash)
        endif()
      endif()

      if(GOOEX_TIME)
        math(EXPR GOOEX_COST "501-${GOOEX_TIME}")
      else()
        set(GOOEX_COST "500")
      endif()

      set_property(TEST ${GOOEX_TEST_NAME} PROPERTY COST ${GOOEX_COST})

    endif()
  endif()
endfunction()

if(ROOT_FOUND)
  add_subdirectory(2d_plot)
  add_subdirectory(addition)
  add_subdirectory(chisquare)
  add_subdirectory(convolution)
  add_subdirectory(dalitz)
  add_subdirectory(pipipi0DPFit)
  add_subdirectory(product)
  add_subdirectory(SigGen)
  add_subdirectory(simpleFit)
  add_subdirectory(zachFit)
  add_subdirectory(DP4)
  add_subdirectory(TDDP4)
endif()

# Root not required
add_subdirectory(exponential)

goofit_add_link(RunAll.py)
