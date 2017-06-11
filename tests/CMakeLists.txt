
goofit_add_library(AllTests_lib AllTests.cu)
set_target_properties(AllTests_lib PROPERTIES FOLDER tests)

goofit_add_executable(AllTests AllTests.cu)
set_target_properties(AllTests PROPERTIES FOLDER tests)

macro(goofit_add_test TESTNAME)
    goofit_add_executable(${TESTNAME} ${TESTNAME}.cu)
    target_link_libraries(${TESTNAME} AllTests_lib)
    target_sources(AllTests PUBLIC ${TESTNAME}.cu)
    add_test(${TESTNAME} ${TESTNAME})
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()

goofit_add_test(BinnedTest)
goofit_add_test(SimpleTest)
goofit_add_test(UnbinnedTest)
#NormalizeTest.cu
#BlindTest.cu)

if(ROOT_FOUND)
    #    list(APPEND TESTS Minuit1Test.cu)
endif()
