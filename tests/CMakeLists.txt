add_subdirectory("${PROJECT_SOURCE_DIR}/extern/googletest" "extern/googletest")

macro(goofit_add_test TESTNAME)
    goofit_add_executable(${TESTNAME} ${TESTNAME}.cu)
    target_link_libraries(${TESTNAME} gtest gmock gtest_main)
    add_test(${TESTNAME} ${TESTNAME})
endmacro()

goofit_add_test(SimpleTest)
goofit_add_test(NormalizeTest)

if(ROOT_FOUND)
    goofit_add_test(Minuit1Test)
endif()
