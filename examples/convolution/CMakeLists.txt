goofit_add_directory()
goofit_add_executable(convolution convolution.cu)
set_target_properties(convolution PROPERTIES FOLDER examples)

add_test(NAME convolution_Example COMMAND convolution)
