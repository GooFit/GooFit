goofit_add_directory()
goofit_add_executable(exponential exponential.cu)
if(ROOT_FOUND)
    goofit_add_executable(exponential1 exponential1.cu)
endif()
goofit_add_executable(exponential2 exponential2.cu)

