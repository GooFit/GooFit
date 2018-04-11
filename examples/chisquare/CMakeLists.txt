goofit_add_directory()
goofit_add_executable(chisquare chisquare.cu)
set_target_properties(chisquare PROPERTIES FOLDER examples)

add_test(NAME chisquare_Example COMMAND chisquare)
