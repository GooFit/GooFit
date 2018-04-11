goofit_add_directory()
goofit_add_executable(dalitz dalitz.cu)
set_target_properties(dalitz PROPERTIES FOLDER examples)

add_test(NAME dalitz_Example COMMAND dalitz)
