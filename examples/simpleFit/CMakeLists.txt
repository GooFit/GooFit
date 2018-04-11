goofit_add_directory()
goofit_add_executable(simpleFit simpleFit.cu)
set_target_properties(simpleFit PROPERTIES FOLDER examples)

add_test(NAME simpleFit_Example COMMAND simpleFit)
