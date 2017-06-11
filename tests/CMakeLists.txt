
goofit_add_executable(AllTests AllTests.cu
    BinnedTest.cu
    SimpleTest.cu
    UnbinnedTest.cu
    )
set_target_properties(AllTests PROPERTIES FOLDER tests)

#NormalizeTest.cu
#BlindTest.cu)

if(ROOT_FOUND)
    #    target_souces(AllTests PUBLIC Minuit1Test.cu)
endif()

add_test(AllTests AllTests)
