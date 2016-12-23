cmake_minimum_required(VERSION 3.4)

project(GOOFIT CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(FATAL_ERROR "You cannot build in the source directory. Please make a build subdirectory.")
endif()

set(DEVICE_LISTING CUDA OMP TBB CPP)
set(HOST_LISTING OMP TBB CPP)

set(GOOFIT_DEVICE CUDA CACHE STRING "The compute device, options are ${DEVICE_LISTING}")
set(GOOFIT_HOST OMP CACHE STRING "The compute device, options are ${HOST_LISTING}")

# Requires CMake 3.3 or greater
if(NOT ${GOOFIT_DEVICE} IN_LIST DEVICE_LISTING)
    message(FATAL_ERROR "You must select a device from ${DEVICE_LISTING}, not ${GOOFIT_DEVICE}")
endif()
if(NOT ${GOOFIT_HOST} IN_LIST HOST_LISTING)
    message(FATAL_ERROR "You must select a host from ${HOST_LISTING}, not ${HOST_DEVICE}")
endif()

add_definitions("-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_${GOOFIT_DEVICE}")
add_definitions("-DMCBOOSTER_BACKEND=${GOOFIT_DEVICE}")
add_definitions("-DTHRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_${GOOFIT_HOST}")

# -Wl,--no-undefined,--no-allow-shlib-undefined")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

if(GOOFIT_DEVICE STREQUAL CUDA)
    find_package(CUDA 6.0 REQUIRED)

    set(CUDA_NVCC_FLAGS; ${CUDA_NVCC_FLAGS}; -std=c++11; ${CMAKE_CXX_FLAGS})
    set(CUDA_SEPARABLE_COMPILATION OFF)
    
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindCudaArch.cmake)
    select_nvcc_arch_flags(NVCC_FLAGS_EXTRA)
    list(APPEND CUDA_NVCC_FLAGS ${NVCC_FLAGS_EXTRA})
else()
    include_directories(${PROJECT_SOURCE_DIR}/include/goofit/fakecuda)
    find_package(Thrust REQUIRED)
    include_directories("${THRUST_INCLUDE_DIRS}")
endif()

if(GOOFIT_DEVICE STREQUAL OMP OR GOOFIT_HOST STREQUAL OMP)
find_package(OpenMP REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
find_package(Threads REQUIRED)
link_libraries(Threads::Threads)
endif()

if(GOOFIT_DEVICE STREQUAL TBB OR GOOFIT_HOST STREQUAL TBB)
find_package(TBB REQUIRED)
include_directories("${TBB_INCLUDE_DIRS}")
endif()


find_package(ROOT 6.08 REQUIRED COMPONENTS Minuit)
include_directories("${ROOT_INCLUDE_DIRS}")
link_libraries(${ROOT_LIBRARIES})


#if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
#elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost -march=native")
#endif()
 
function(goofit_add_library GNAME)
    if(GOOFIT_DEVICE STREQUAL CUDA)
        cuda_add_library(${GNAME} STATIC ${ARGN})
    else()
        foreach(N ${ARGN})
            SET_SOURCE_FILES_PROPERTIES(${N} PROPERTIES LANGUAGE CXX)
        endforeach()
        add_library(${GNAME} STATIC ${ARGN})
        target_compile_options(${GNAME} PUBLIC -x c++)
    endif()
endfunction()

function(goofit_add_executable NAMEEXE)
    if(GOOFIT_DEVICE STREQUAL CUDA)
        cuda_add_executable(${NAMEEXE} ${ARGN})
        link_libraries(${CUDA_LIBRARIES})
    else()
        foreach(N ${ARGN})
            SET_SOURCE_FILES_PROPERTIES(${N} PROPERTIES LANGUAGE CXX)
        endforeach()
	    add_executable(${NAMEEXE} ${ARGN})
    endif()
    target_link_libraries(${NAMEEXE} goofit_lib)
    target_compile_options(${NAMEEXE} PUBLIC -x c++)
endfunction()

add_definitions("-DTARGET_SM35")
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/MCBooster)

add_subdirectory(src)

add_library(goofit_lib INTERFACE)
target_link_libraries(goofit_lib INTERFACE goofit_base PDFs rootstuff)

add_subdirectory(examples)

export(PACKAGE GooFit)
