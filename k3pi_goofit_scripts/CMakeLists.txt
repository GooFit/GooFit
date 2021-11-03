cmake_minimum_required(VERSION 3.9...3.15)

if(CMAKE_PROJECT_NAME STREQUAL "GOOFIT")
    goofit_add_package(K3PI)
else()
    project(K3Pi LANGUAGES CXX)
    add_subdirectory(goofit)
    goofit_external_package()
endif()

goofit_add_executable(amp4body_cpp cpp_amp4body_TD_example.cu)
goofit_add_executable(amp4body_cpp_py cpp_python_amp4body_td_example.cu)
target_link_libraries(amp4body_cpp_py PRIVATE pybind11::embed)
