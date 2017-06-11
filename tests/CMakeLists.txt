
goofit_add_library(AllTests AllTests.cu)
set_target_properties(AllTests PROPERTIES FOLDER tests)

macro(goofit_add_test TESTNAME)
    goofit_add_executable(${TESTNAME} ${TESTNAME}.cu)
    target_link_libraries(${TESTNAME} AllTests)
    add_test(${TESTNAME} ${TESTNAME})
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()

goofit_add_test(BinnedTest)
goofit_add_test(SimpleTest)
#UnbinnedTest.cu
#NormalizeTest.cu
#BlindTest.cu)

if(ROOT_FOUND)
    #    list(APPEND TESTS Minuit1Test.cu)
endif()
