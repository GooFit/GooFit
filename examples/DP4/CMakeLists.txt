goofit_add_directory()
goofit_add_executable(DP4 DP4.cu)
set_target_properties(DP4 PROPERTIES FOLDER examples)

add_test(NAME DP4_Example COMMAND DP4)
