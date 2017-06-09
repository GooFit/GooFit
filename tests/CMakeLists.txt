add_subdirectory("${PROJECT_SOURCE_DIR}/extern/googletest" "extern/googletest")

mark_as_advanced(
    BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
    gmock_build_tests gtest_build_samples gtest_build_tests
    gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)

set_target_properties(gtest PROPERTIES FOLDER extern)
set_target_properties(gtest_main PROPERTIES FOLDER extern)
set_target_properties(gmock PROPERTIES FOLDER extern)
set_target_properties(gmock_main PROPERTIES FOLDER extern)

macro(goofit_add_test TESTNAME)
    goofit_add_executable(${TESTNAME} ${TESTNAME}.cu)
    target_link_libraries(${TESTNAME} gtest gmock gtest_main)
    add_test(${TESTNAME} ${TESTNAME})
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()

goofit_add_test(SimpleTest)
goofit_add_test(BinnedTest)
goofit_add_test(UnbinnedTest)
goofit_add_test(NormalizeTest)
goofit_add_test(BlindTest)

if(ROOT_FOUND)
    goofit_add_test(Minuit1Test)
endif()
