goofit_add_directory()
goofit_add_executable(addition addition.cu)
set_target_properties(addition PROPERTIES FOLDER examples)

add_test(NAME addition_Example COMMAND addition)
