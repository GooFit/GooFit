goofit_add_directory()
goofit_add_executable(2d_plot 2d_plot.cu)
set_target_properties(2d_plot PROPERTIES FOLDER examples)

add_test(NAME 2d_plot_Example COMMAND 2d_plot)
