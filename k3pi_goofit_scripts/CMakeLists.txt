cmake_minimum_required(VERSION 3.9...3.15)

goofit_add_executable(amp4body_cpp cpp_amp4body_TD_example.cu)
goofit_add_executable(amp4body_cpp_py cpp_python_amp4body_td_example.cu)
target_link_libraries(amp4body_cpp_py PRIVATE pybind11::embed)
