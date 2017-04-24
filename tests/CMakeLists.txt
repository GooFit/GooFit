add_subdirectory(${PROJECT_SOURCE_DIR}/extern/googletest)

macro(add_gtest TESTNAME)
    target_link_libraries(${TESTNAME} PUBLIC gtest gmock gtest_main)
    add_test(${TESTNAME} ${TESTNAME})
endmacro()

goofit_add_executable(SimpleTest SimpleTest.cu)
add_gtest(SimpleTest)
