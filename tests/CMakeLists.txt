
if(ROOT_FOUND)
    set(ROOT_LISTING Minuit1Test.cu)
endif()

goofit_add_executable(AllTests AllTests.cu
    BinnedTest.cu
    SimpleTest.cu
    UnbinnedTest.cu
    BlindTest.cu
    NormalizeTest.cu
    ${ROOT_LISTING}
    )

set_target_properties(AllTests PROPERTIES FOLDER tests)

add_test(AllTests AllTests)
