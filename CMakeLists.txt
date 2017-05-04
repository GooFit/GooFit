cmake_minimum_required(VERSION 3.4)

project(GOOFIT
    VERSION 2.0.0
    LANGUAGES CXX)

if(NOT DEFINED CMAKE_SUPPRESS_DEVELOPER_WARNINGS)
    set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE INTERNAL "No dev warnings")
endif()

### Require out-of-source builds except in main directory
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(NOT PROJECT_BINARY_DIR STREQUAL PROJECT_SOURCE_DIR AND EXISTS ${LOC_PATH})
    message(FATAL_ERROR "You cannot build in a sub-source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory.")
endif()

# Update submodules as needed
option(GOOFIT_SUBMODULE "Check submodules during build" ON)
if(GOOFIT_SUBMODULE)
    message(STATUS "Submodule update")
    execute_process(COMMAND git submodule update --init --recursive
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                    RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
        message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
endif()

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Add support for better CUDA behaviour in CMake < 3.7
if(CMAKE_VERSION VERSION_LESS 3.7)
    set(CMAKE_MODULE_PATH  "${PROJECT_SOURCE_DIR}/extern/support_old_cmake" ${CMAKE_MODULE_PATH})
endif()

# Add Sanatizers
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/extern/sanitizers/cmake" ${CMAKE_MODULE_PATH})
find_package(Sanitizers)


# Output the current GooFit version
configure_file (
    "${PROJECT_SOURCE_DIR}/include/goofit/Version.h.in"
    "${PROJECT_BINARY_DIR}/include/goofit/Version.h"
)


set(GOOFIT_CUDA_OR_GROUPSIZE "128" CACHE STRING "Overrides the default group distribution for Thrust's transform_reduce")
set(GOOFIT_CUDA_OR_GRAINSIZE "7" CACHE STRING "Overrides the default grain size for Thrust's transform_reduce")

configure_file (
    "${PROJECT_SOURCE_DIR}/include/goofit/ThrustOverrideConfig.h.in"
    "${PROJECT_BINARY_DIR}/include/goofit/ThrustOverrideConfig.h"
)

include_directories("${PROJECT_BINARY_DIR}/include")

### C++ settings ###
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_CXX_VISIBILITY_PRESET hidden)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
endif()

# Printout time taken to compile
option(GOOFIT_TIME_COMPILE "Print time to compile during compilation" OFF)
if(GOOFIT_TIME_COMPILE)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "time")
endif()

# Code coverage
if(CMAKE_BUILD_TYPE STREQUAL Coverage)
    include(CodeCoverage)
    setup_target_for_coverage(GooFit_coverage ctest coverage)
endif()

### Options ###
set(DEVICE_LISTING CUDA OMP CPP TBB Auto)
set(HOST_LISTING OMP CPP TBB Auto)
mark_as_advanced(DEVICE_LISTING HOST_LISTING)

set(GOOFIT_DEVICE Auto CACHE STRING "The compute device, options are ${DEVICE_LISTING}")
set(GOOFIT_HOST Auto CACHE STRING "The compute device, options are ${HOST_LISTING}")

# Requires CMake 3.3 or greater
if(NOT ${GOOFIT_DEVICE} IN_LIST DEVICE_LISTING)
    message(FATAL_ERROR "You must select a device from ${DEVICE_LISTING}, not ${GOOFIT_DEVICE}")
endif()
if(NOT ${GOOFIT_HOST} IN_LIST HOST_LISTING)
    message(FATAL_ERROR "You must select a host from ${HOST_LISTING}, not ${HOST_DEVICE}")
endif()

# Auto device selection
if(GOOFIT_DEVICE STREQUAL Auto)
    find_package(CUDA 7.0)
    if(CUDA_FOUND)
        set(GOOFIT_DEVICE CUDA)
    else()
        find_package(OpenMP QUIET)
        if(OpenMP_FOUND)
            set(GOOFIT_DEVICE OMP)
        else()
            set(GOOFIT_DEVICE CPP)
        endif()
    endif()
    message(STATUS "Auto device selection: ${GOOFIT_DEVICE}")
endif()

# Auto host selection based on device
if(GOOFIT_HOST STREQUAL Auto)
    if(GOOFIT_DEVICE STREQUAL OMP)
        set(GOOFIT_HOST OMP)
    elseif(GOOFIT_DEVICE STREQUAL TBB)
        set(GOOFIT_HOST TBB)
    else()
        set(GOOFIT_HOST CPP)
    endif()
endif()

# Checks for invalid combinations
if(${GOOFIT_DEVICE} STREQUAL TBB AND ${GOOFIT_HOST} STREQUAL OMP)
    message(FATAL_ERROR "You must set TBB as both host and device (OMP will still be required)")
endif()
if(${GOOFIT_DEVICE} STREQUAL OMP AND ${GOOFIT_HOST} STREQUAL TBB)
    message(FATAL_ERROR "TBB cannot be a host backend for OMP")
endif()
if(${GOOFIT_DEVICE} STREQUAL CUDA AND ${GOOFIT_HOST} STREQUAL OMP)
    message(FATAL_ERROR "OMP cannot be a host backend for CUDA")
endif()

# Only default to off for CUDA
if(GOOFIT_DEVICE STREQUAL CUDA)
    option(GOOFIT_SEPARATE_COMP "Enable separate compilation of PDFs" OFF)
else()
    option(GOOFIT_SEPARATE_COMP "Enable separate compilation of PDFs" ON)
endif()

if(GOOFIT_SEPARATE_COMP)
    add_definitions("-DSEPARABLE")
endif()

add_definitions("-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_${GOOFIT_DEVICE}")
add_definitions("-DMCBOOSTER_BACKEND=${GOOFIT_DEVICE}")
add_definitions("-DTHRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_${GOOFIT_HOST}")

option(GOOFIT_DEBUG "Print debugging messages" OFF)
option(GOOFIT_TRACE "Print messages to trace the behavior of GooFit" OFF)
if(GOOFIT_DEBUG)
    add_definitions("-DGOOFIT_DEBUG_FLAG=1")
endif()
if(GOOFIT_TRACE)
    add_definitions("-DGOOFIT_TRACE_FLAG=1")
endif()


set(GOOFIT_ARCH Auto CACHE STRING "The GPU Archetecture, can be Auto, All, Common, a number, or a name")

option(GOOFIT_MPI "Turn on MPI for goofit" OFF)

if(GOOFIT_MPI)
    find_package(MPI REQUIRED)

    # Added globally
    add_definitions("-DGOOFIT_MPI")
    add_compile_options(${MPI_CXX_COMPILE_FLAGS})
    include_directories(${MPI_CXX_INCLUDE_PATH})
    list(APPEND CMAKE_EXE_LINKER_FLAGS ${MPI_CXX_LINK_FLAGS})
    link_libraries(${MPI_CXX_LIBRARIES})

    message(STATUS "MPI found. Use the following to run your program")
    message(STATUS "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} PROCS ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
    message(STATUS "where PROCS is the number of processors on which to execute the program, EXECUTABLE is the MPI program, and ARGS are the arguments to pass to the MPI program.")
endif()

if(GOOFIT_DEVICE STREQUAL CUDA)
    if(NOT CUDA_FOUND)
        find_package(CUDA 7.0 REQUIRED)
    endif()

    list(APPEND CUDA_NVCC_FLAGS -std=c++11)
    list(APPEND CUDA_NVCC_FLAGS ${MPI_CXX_LINK_FLAGS})

    if(GOOFIT_SEPARATE_COMP)
        set(CUDA_SEPARABLE_COMPILATION ON)
    else()
        set(CUDA_SEPARABLE_COMPILATION OFF)
    endif()

    cuda_select_nvcc_arch_flags(ARCH_FLAGS ${GOOFIT_ARCH})
    list(APPEND CUDA_NVCC_FLAGS ${ARCH_FLAGS})
    message(STATUS "Compiling for GPU arch: ${ARCH_FLAGS}")

endif()

set(GOOFIT_CFLAGS "$<$<CONFIG:Release>:-march=native>$<$<CONFIG:RelWithDebInfo>:-march=native>" CACHE STRING "compiler flags")
add_compile_options("${GOOFIT_CFLAGS}")

set(THRUST_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/extern/thrust")
find_package(Thrust 1.8 REQUIRED)
include_directories(SYSTEM "${THRUST_INCLUDE_DIRS}")

if(GOOFIT_DEVICE STREQUAL OMP OR GOOFIT_HOST STREQUAL OMP OR GOOFIT_DEVICE STREQUAL TBB OR GOOFIT_HOST STREQUAL TBB)
    find_package(OpenMP REQUIRED)
    add_compile_options(${OpenMP_CXX_FLAGS})
    find_package(Threads REQUIRED)
    link_libraries(Threads::Threads)
    link_libraries(${OpenMP_CXX_FLAGS})
endif()

if(GOOFIT_DEVICE STREQUAL TBB OR GOOFIT_HOST STREQUAL TBB)
    find_package(TBB COMPONENTS tbbmalloc tbbmalloc_proxy tbb_preview)
    #set(CMAKE_CXX_FLAGS -xMIC-AVX512)
    #set(CMAKE_CXX_FLAGS "-xHOST -qopenmp")
    include_directories(SYSTEM "${TBB_INCLUDE_DIRS}")
    link_libraries(${TBB_LIBRARIES})
endif()

# Include directories are not picked up by FindCUDA
find_package(ROOT 6 COMPONENTS Minuit Minuit2)
if(NOT ROOT_FOUND)
    add_subdirectory("extern/Minuit2")
endif()

# Adding simple libraries
add_subdirectory("extern/CLI11")
add_subdirectory("extern/FeatureDetector")
add_subdirectory("extern/fmt")

add_library(rang INTERFACE)
target_include_directories(rang INTERFACE "${PROJECT_SOURCE_DIR}/extern/rang/include")

function(goofit_add_library GNAME)
    if(GOOFIT_DEVICE STREQUAL CUDA)
        foreach(N ${ARGN})
            get_filename_component(NEXT ${N} EXT)
            if(NEXT STREQUAL ".cc" OR NEXT STREQUAL ".cpp" OR NEXT STREQUAL ".cxx")
                set_source_files_properties(${N} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
            endif()
        endforeach()
        cuda_add_library(${GNAME} STATIC ${ARGN} OPTIONS "-Xcompiler=-fPIC")
    else()
        foreach(N ${ARGN})
            get_filename_component(NEXT ${N} EXT)
            if(NEXT STREQUAL ".cu")
                set_source_files_properties(${N} PROPERTIES LANGUAGE CXX)
            endif()
        endforeach()
        add_library(${GNAME} STATIC ${ARGN})
        add_sanitizers(${GNAME})
        target_compile_options(${GNAME} PUBLIC -x c++)
    endif()
    target_link_libraries(${GNAME} rang CLI11 fmt)

    if(ROOT_FOUND)
        target_link_libraries(${GNAME} ROOT::ROOT)
    else()
        target_link_libraries(${GNAME} Minuit2)
    endif()
endfunction()

function(goofit_add_executable NAMEEXE)
    if(GOOFIT_DEVICE STREQUAL CUDA)
        cuda_add_executable(${NAMEEXE} ${ARGN})
        link_libraries(${CUDA_LIBRARIES})
    else()
        foreach(N ${ARGN})
            get_filename_component(NEXT ${N} EXT)
            if(NEXT STREQUAL ".cu")
                set_source_files_properties(${N} PROPERTIES LANGUAGE CXX)
            endif()
        endforeach()
        add_executable(${NAMEEXE} ${ARGN})
        add_sanitizers(${NAMEEXE})
    endif()
    target_link_libraries(${NAMEEXE} goofit_lib)
    target_compile_options(${NAMEEXE} PUBLIC -x c++)
endfunction()

function(goofit_add_link)
    foreach(NAMELINK ${ARGN})

        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${NAMELINK}")
            add_custom_command(
                OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${NAMELINK}"
                COMMAND ${CMAKE_COMMAND} -E create_symlink
                "${CMAKE_CURRENT_SOURCE_DIR}/${NAMELINK}"
                "${CMAKE_CURRENT_BINARY_DIR}/${NAMELINK}"
                )

            get_filename_component(a_dir ${CMAKE_CURRENT_SOURCE_DIR} NAME)

            add_custom_target(
                "${a_dir}_${NAMELINK}_link" ALL DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${NAMELINK}"
                )
        else()
            message(STATUS "Warning: ${NAMELINK} does not exist. Not making link. Build then retry?")
        endif()
    endforeach()
endfunction()


function(goofit_add_directory)
    file(GLOB directory_listing RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *)
    set(skip_files CMakeLists.txt Makefile makefile .gitignore .git)

    foreach(NAMELINK ${directory_listing})
        if(NOT ${NAMELINK} IN_LIST skip_files)
            add_custom_command(
                OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${NAMELINK}"
                COMMAND ${CMAKE_COMMAND} -E create_symlink
                "${CMAKE_CURRENT_SOURCE_DIR}/${NAMELINK}"
                "${CMAKE_CURRENT_BINARY_DIR}/${NAMELINK}"
                )

            list(APPEND dir_listing "${CMAKE_CURRENT_BINARY_DIR}/${NAMELINK}")

        endif()
    endforeach()

    get_filename_component(a_dir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    add_custom_target("${a_dir}_directory" ALL DEPENDS ${dir_listing})
endfunction()


# Adding global property to collect extra info for compilation
set_property(GLOBAL PROPERTY GOOFIT_EXTRA_PDFS "")
macro(goofit_add_pdf)
    foreach(N ${ARGN})
        get_filename_component(N ${N} ABSOLUTE)
        message(STATUS "Adding PDF: ${N}")
        get_property(goo_extra GLOBAL PROPERTY GOOFIT_EXTRA_PDFS)
        set_property(GLOBAL PROPERTY GOOFIT_EXTRA_PDFS ${goo_extra} ${N})
    endforeach()
endmacro()

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/extern/MCBooster)
include_directories(${PROJECT_SOURCE_DIR}/extern/generics)

add_subdirectory(src)

add_library(goofit_lib INTERFACE)
target_link_libraries(goofit_lib INTERFACE goofit_base PDFs rang CLI11 FeatureDetector)

if(ROOT_FOUND)
    target_link_libraries(goofit_lib INTERFACE ROOT::ROOT ROOT::Minuit ROOT::Minuit2)
else()
    target_link_libraries(goofit_lib INTERFACE Minuit2)
endif()


option(GOOFIT_EXAMPLES "Build the example programs" ON)
if(GOOFIT_EXAMPLES)
    add_subdirectory(examples)
endif()

if(EXISTS work)
    add_subdirectory(work)
endif()

option(GOOFIT_PACKAGES "Build any goofit* packages found" ON)
if(GOOFIT_PACKAGES)
    file(GLOB list_of_packages RELATIVE ${PROJECT_SOURCE_DIR} goofit*)
    foreach(d ${list_of_packages})
        add_subdirectory(${d})
    endforeach()
endif()


option(GOOFIT_PYTHON "Python bindings for goofit" OFF)
if(GOOFIT_PYTHON)
    add_subdirectory("extern/pybind11")
    message(STATUS "Found Python at ${PYTHON_LIBRARIES}")
    add_subdirectory(python)
endif()

option(GOOFIT_TESTS "Build the goofit tests" OFF)
if(GOOFIT_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

