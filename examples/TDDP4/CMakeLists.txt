goofit_add_directory()
goofit_add_executable(TDDP4 TDDP4.cu)
set_target_properties(TDDP4 PROPERTIES FOLDER examples)

add_test(NAME TDDP4_Example COMMAND TDDP4 -t 5)
