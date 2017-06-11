include(cmake/ParseAndAddCatchTests.cmake)

set(TESTS
BinnedTest.cu)
#SimpleTest.cu
#UnbinnedTest.cu
#NormalizeTest.cu
#BlindTest.cu)

if(ROOT_FOUND)
    #    list(APPEND TESTS Minuit1Test.cu)
endif()

goofit_add_executable(AllTests AllTests.cu ${TESTS})
ParseAndAddCatchTests(AllTests)
set_target_properties(AllTests PROPERTIES FOLDER tests)
