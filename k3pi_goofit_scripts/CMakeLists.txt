cmake_minimum_required(VERSION 3.9...3.15)

if(CMAKE_PROJECT_NAME STREQUAL "GOOFIT")
    goofit_add_package(K3PI)
else()
    project(K3Pi LANGUAGES CXX)
    add_subdirectory(goofit)
    goofit_external_package()
endif()

goofit_add_executable(amp4body_cpp cpp_amp4body_TD_example.cu)

if(GOOFIT_KMATRIX)
    message(STATUS "K3pi: Building kMatrix based formalism")
    # Not needed in GooFit 2.2
    set(GOOFIT_MAXPAR 5900 CACHE STRING "Required by K3pi" FORCE)
    set(GOOFIT_MAXFUNC 300 CACHE STRING "Required by K3pi" FORCE)
    set_property(SOURCE ampgen/DtoKpipipi_v2.cu PROPERTY HEADER_FILE_ONLY ON)
    goofit_add_executable(K3pi_v2 K3pi_v2.cu ampgen/DtoKpipipi_v2.cu)
else()
    message(STATUS "K3pi: NOT building kMatrix based formalism, use -DGOOFIT_KMATRIX=ON to build (GPU not supported yet)")
endif()