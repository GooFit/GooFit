cmake_minimum_required(VERSION 3.4)

project(GOOFIT
    VERSION 2.1.0
    LANGUAGES CXX)

set(GOOFIT_TAG "dev")
#set(GOOFIT_TAG "alpha")
#set(GOOFIT_TAG "beta")
#set(GOOFIT_TAG "release")

if(NOT DEFINED CMAKE_SUPPRESS_DEVELOPER_WARNINGS AND NOT GOOFIT_TAG STREQUAL "dev")
    set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE INTERNAL "No dev warnings")
endif()

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS ${LOC_PATH})
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

# Allow IDE's to group targets into folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Get the git command
find_package(Git QUIET)

if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GOOFIT_SUBMODULE "Check submodules during build" ON)
    if(GOOFIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()

    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
                    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                    OUTPUT_VARIABLE GOOFIT_GIT_VERSION
                    ERROR_QUIET
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

else()
    set(GOOFIT_GIT_VERSION "unknown")
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/cmake/FindThrust.cmake"
OR NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/CLI11/CMakeLists.txt"
OR NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/Eigen/CMakeLists.txt"
OR NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/FeatureDetector/CMakeLists.txt"
OR NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/MCBooster/CMakeLists.txt"
OR NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/thrust/README.md")
    message(FATAL_ERROR "The submodules were not downloaded! GOOFIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Add support for better CUDA behaviour in CMake < 3.7
if(CMAKE_VERSION VERSION_LESS 3.7)
    set(CMAKE_MODULE_PATH  "${PROJECT_SOURCE_DIR}/extern/support_old_cmake" ${CMAKE_MODULE_PATH})
endif()

if(CMAKE_VERSION VERSION_GREATER 3.6)
# Add clang-tidy if available
option(GOOFIT_TIDY_FIX "Perform fixes for Clang-Tidy - changes source inplace" OFF)
    find_program(
        CLANG_TIDY_EXE
        NAMES "clang-tidy"
        DOC "Path to clang-tidy executable"
    )

    if(CLANG_TIDY_EXE)
        if(GOOFIT_TIDY_FIX)
            set(DO_CLANG_TIDY "${CLANG_TIDY_EXE}" "-fix")
        else()
            set(DO_CLANG_TIDY "${CLANG_TIDY_EXE}")
        endif()
    endif()
endif()

# Add Sanatizers
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/extern/sanitizers/cmake" ${CMAKE_MODULE_PATH})
find_package(Sanitizers)

set(GOOFIT_MAXPAR "1800" CACHE STRING "The number of parameters to statically support, can be increased but should not be too large.")

set(GOOFIT_CUDA_OR_GROUPSIZE "128" CACHE STRING "Overrides the default group distribution for Thrust's transform_reduce")
set(GOOFIT_CUDA_OR_GRAINSIZE "7" CACHE STRING "Overrides the default grain size for Thrust's transform_reduce")

configure_file (
    "${PROJECT_SOURCE_DIR}/include/goofit/detail/ThrustOverrideConfig.h.in"
    "${PROJECT_BINARY_DIR}/include/goofit/detail/ThrustOverrideConfig.h"
)

include_directories("${PROJECT_BINARY_DIR}/include")

### C++ settings ###
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_CXX_VISIBILITY_PRESET hidden)

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unknown-pragmas -Wno-long-long -Wno-attributes -Wno-sign-compare -Wno-unused-parameter")
    # Helpful but irritating: -Wzero-as-null-pointer-constant -Wsuggest-override
    # no-sign-compare can be removed, but will take some work to clean up
    # Same is true for no-unused-parameter
    # if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.0)
endif()

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

option(GOOFIT_SEPARATE_COMP "Enable separate compilation of PDFs" ON)

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

# Adding backtrace (optional)
# Some systems need execinfo explicitly linked
# Standard CMake module
find_package(Backtrace)
add_library(backtrace INTERFACE)
if(Backtrace_FOUND)
    # Assuming no extra flags
    target_include_directories(backtrace INTERFACE ${Backtrace_INCLUDE_DIR})
    target_link_libraries(backtrace INTERFACE ${Backtrace_LIBRARIES})
endif()
configure_file(
    "${PROJECT_SOURCE_DIR}/include/goofit/detail/Backtrace.h.in"
    "${PROJECT_BINARY_DIR}/include/goofit/detail/Backtrace.h"
)

set(GOOFIT_ARCH Auto CACHE STRING "The GPU Archetecture, can be Auto, All, Common, a number, or a name")

option(GOOFIT_MPI "Turn on MPI for goofit" OFF)


if(GOOFIT_MPI)
    find_package(MPI REQUIRED)

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

    if(GOOFIT_MPI)
        list(APPEND CUDA_NVCC_FLAGS ${MPI_CXX_COMPILE_FLAGS})
    endif()

    list(APPEND CUDA_NVCC_FLAGS -std=c++11)

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

# This library provides an interface to openMP, etc for all libraries.
add_library(goofit_mt INTERFACE)

if(GOOFIT_DEVICE STREQUAL OMP OR GOOFIT_HOST STREQUAL OMP OR GOOFIT_DEVICE STREQUAL TBB OR GOOFIT_HOST STREQUAL TBB)
    find_package(OpenMP REQUIRED)
    target_compile_options(goofit_mt INTERFACE ${OpenMP_CXX_FLAGS})
    find_package(Threads REQUIRED)
    target_link_libraries(goofit_mt INTERFACE Threads::Threads)
    target_link_libraries(goofit_mt INTERFACE ${OpenMP_CXX_FLAGS})
endif()

if(GOOFIT_DEVICE STREQUAL TBB OR GOOFIT_HOST STREQUAL TBB)
    find_package(TBB COMPONENTS tbbmalloc tbbmalloc_proxy tbb_preview)
    target_include_directories(goofit_mt INTERFACE SYSTEM "${TBB_INCLUDE_DIRS}")
    target_link_libraries(goofit_mt INTERFACE ${TBB_LIBRARIES})
endif()

# Include directories are not picked up by FindCUDA
find_package(ROOT 6 COMPONENTS Minuit)
if(ROOT_FOUND)
    add_definitions("-DROOT_FOUND")
endif()
# If either ROOT not found or Minuit2 not included
if(NOT ROOT_Minuit2_LIBRARY)
    add_subdirectory("extern/Minuit2")
    set_target_properties(Minuit2 PROPERTIES FOLDER extern)
    set_target_properties(Math PROPERTIES FOLDER extern)
endif()

# Adding simple libraries
add_subdirectory("extern/CLI11")
mark_as_advanced(CLI_EXAMPLES CLI_SINGLE_FILE CLI_SINGLE_FILE_TESTS CLI_TESTING)
add_subdirectory("extern/FeatureDetector")
set_target_properties(FeatureDetector PROPERTIES FOLDER extern)

## Format
add_subdirectory("extern/fmt")
set_target_properties(fmt PROPERTIES FOLDER extern)
mark_as_advanced(FMT_CPPFORMAT FMT_DOC FMT_INSTALL FMT_PEDANTIC FMT_TEST FMT_USE_CPP11)


add_library(rang INTERFACE)
target_include_directories(rang INTERFACE "${PROJECT_SOURCE_DIR}/extern/rang/include")


add_library(Eigen INTERFACE)
target_include_directories(Eigen INTERFACE "${PROJECT_SOURCE_DIR}/extern/Eigen/Eigen")

# Output the current GooFit version
configure_file (
    "${PROJECT_SOURCE_DIR}/include/goofit/Version.h.in"
    "${PROJECT_BINARY_DIR}/include/goofit/Version.h"
)


function(goofit_add_library GNAME)
    set(options STATIC SHARED MODULE)
    set(oneValueArgs "")
    set(multiValueArgs "")
    cmake_parse_arguments(
        GOOFIT_LIB
        "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN}
        )
    if(GOOFIT_LIB_SHARED)
        set(GOOFIT_LIB_TYPE SHARED)
    elseif(GOOFIT_LIB_MODULE)
        set(GOOFIT_LIB_TYPE MODULE)
    else()
        set(GOOFIT_LIB_TYPE STATIC)
    endif()

    if(GOOFIT_DEVICE STREQUAL CUDA)
        foreach(N ${GOOFIT_LIB_UNPARSED_ARGUMENTS})
            get_filename_component(NEXT ${N} EXT)
            if(NEXT STREQUAL ".cc" OR NEXT STREQUAL ".cpp" OR NEXT STREQUAL ".cxx")
                set_source_files_properties(${N} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
            endif()
        endforeach()
        cuda_add_library(${GNAME} ${GOOFIT_LIB_TYPE} ${GOOFIT_LIB_UNPARSED_ARGUMENTS} OPTIONS "-Xcompiler=-fPIC")
    else()
        foreach(N ${GOOFIT_LIB_UNPARSED_ARGUMENTS})
            get_filename_component(NEXT ${N} EXT)
            if(NEXT STREQUAL ".cu")
                set_source_files_properties(${N} PROPERTIES LANGUAGE CXX)
            endif()
        endforeach()
        add_library(${GNAME} ${GOOFIT_LIB_TYPE} ${GOOFIT_LIB_UNPARSED_ARGUMENTS})
        add_sanitizers(${GNAME})
        target_compile_options(${GNAME} PUBLIC -x c++)
    endif()
    target_link_libraries(${GNAME} rang CLI11 fmt goofit_mt)
    set_target_properties(${GNAME} PROPERTIES FOLDER core)

    set_target_properties(${GNAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION 1)

    if(ROOT_FOUND)
        target_link_libraries(${GNAME} ROOT::ROOT)
    else()
        target_link_libraries(${GNAME} Minuit2)
    endif()

    source_group("Source Files" REGULAR_EXPRESSION ".*\\.c[uc]")

    if(CLANG_TIDY_EXE)
        set_target_properties(
            ${GNAME} PROPERTIES
            CXX_CLANG_TIDY "${DO_CLANG_TIDY}"
        )
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
    set_target_properties(${NAMEEXE} PROPERTIES FOLDER projects)
    source_group("Source Files" REGULAR_EXPRESSION ".*\\.c[uc]")

    set_target_properties(${GNAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION 1)

    if(CLANG_TIDY_EXE)
        set_target_properties(
            ${GNAME} PROPERTIES
            CXX_CLANG_TIDY "${DO_CLANG_TIDY}"
        )
    endif()
endfunction()

function(goofit_add_link)
    if(MSVC) # Not officially supported, but needed to even configure on Windows
        set(COMM copy)
    else()
        set(COMM create_symlink)
    endif()

    foreach(NAMELINK ${ARGN})
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${NAMELINK}")

            execute_process(
                COMMAND ${CMAKE_COMMAND} -E ${COMM}
                "${CMAKE_CURRENT_SOURCE_DIR}/${NAMELINK}"
                "${CMAKE_CURRENT_BINARY_DIR}/${NAMELINK}"
                WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
                )

        else()
            message(WARNING "${NAMELINK} does not exist. Not making link.")
        endif()
    endforeach()
endfunction()


function(goofit_add_directory)
    file(GLOB directory_listing RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *)
    set(skip_files CMakeLists.txt CMakeFiles Makefile makefile .gitignore .git)

    foreach(NAMELINK ${directory_listing})
        if(NOT ${NAMELINK} IN_LIST skip_files)
            goofit_add_link(${NAMELINK})
        endif()
    endforeach()
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

option(GOOFIT_PACKAGES "Build any goofit_* packages found" ON)
if(GOOFIT_PACKAGES)
    file(GLOB list_of_packages RELATIVE ${PROJECT_SOURCE_DIR} goofit_*)
    foreach(d ${list_of_packages})
        add_subdirectory(${d})
    endforeach()
endif()

# This allows GOOFIT_PYTHON to be automatic
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/extern/pybind11/tools" ${CMAKE_MODULE_PATH})
find_package(PythonLibsNew)
if(PYTHONLIBS_FOUND)
    option(GOOFIT_PYTHON "Python bindings for goofit" ON)
else()
    option(GOOFIT_PYTHON "Python bindings for goofit" OFF)
endif()

if(GOOFIT_PYTHON)
    add_subdirectory("extern/pybind11")
    message(STATUS "Found Python at ${PYTHON_LIBRARIES}, building bindings")
    add_subdirectory(python)
endif()

option(GOOFIT_TESTS "Build the goofit tests" ON)
if(GOOFIT_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

