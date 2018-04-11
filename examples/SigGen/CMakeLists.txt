goofit_add_directory()
goofit_add_executable(SigGen SigGen.cu)
set_target_properties(SigGen PROPERTIES FOLDER examples)

add_test(NAME SigGen_Example COMMAND SigGen)
