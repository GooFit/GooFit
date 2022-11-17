goofit_add_directory()

goofit_add_library(ToyModel ToyModel.cu ToyModel.h)

goofit_add_executable(GenerateToy GenerateToy.cu)
target_link_libraries(GenerateToy PUBLIC ToyModel)

goofit_add_executable(FitToy FitToy.cu)
target_link_libraries(FitToy PUBLIC ToyModel)
