add_definitions("-DHAVE_ROOT")
goofit_add_executable(SimpleFit simpleFit.cu)
link_libraries(goofit_lib)

