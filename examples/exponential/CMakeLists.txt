goofit_add_directory()
goofit_add_executable(exponential exponential.cu)
set_target_properties(exponential PROPERTIES FOLDER examples)
add_test(NAME Exponential_Example COMMAND exponential)

if(ROOT_FOUND)
    goofit_add_executable(exponential1 exponential1.cu)
    set_target_properties(exponential1 PROPERTIES FOLDER examples)
    add_test(NAME Exponential1_Example COMMAND exponential1)
endif()

goofit_add_executable(exponential2 exponential2.cu)
set_target_properties(exponential2 PROPERTIES FOLDER examples)
add_test(NAME Exponential2_Example COMMAND exponential2)

