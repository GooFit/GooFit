cmake_minimum_required(VERSION 3.16...3.23)

project(
  GOOFIT
  VERSION 2.3.0
  LANGUAGES CXX)

#set(GOOFIT_TAG "dev")
#set(GOOFIT_TAG "alpha")
#set(GOOFIT_TAG "beta")
set(GOOFIT_TAG "release")

set(CLI11_PRECOMPILED ON)

include(FeatureSummary)
include(CMakePrintHelpers)

macro(FEATURE_OPTION NAME)
  get_property(
    MSG
    CACHE ${NAME}
    PROPERTY HELPSTRING)
  add_feature_info(${NAME} ${NAME} "${MSG}")
endmacro()

### Require out-of-source builds
file(TO_CMAKE_PATH "${GOOFIT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
  message(
    FATAL_ERROR
      "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles."
  )
endif()

# Allow IDE's to group targets into folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Check to see if this is the master project
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(CUR_PROJ ON)
else()
  set(CUR_PROJ OFF)
endif()

# Get the git command
find_package(Git QUIET)

set(GIT_DIR "${GOOFIT_SOURCE_DIR}/.git")
if(GIT_FOUND AND EXISTS "${GIT_DIR}")
  # Update submodules as needed
  option(GOOFIT_SUBMODULE "Check submodules during build" OFF)
  if(GOOFIT_SUBMODULE)
    message(STATUS "Submodule update")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} --git-dir="${GIT_DIR}" submodule update --init
      WORKING_DIRECTORY ${GOOFIT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(
        FATAL_ERROR
          "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules"
      )
    endif()
  endif()

  execute_process(
    COMMAND ${GIT_EXECUTABLE} --git-dir="${GIT_DIR}" rev-parse --short HEAD
    WORKING_DIRECTORY "${GOOFIT_SOURCE_DIR}"
    OUTPUT_VARIABLE GOOFIT_GIT_VERSION
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

else()
  set(GOOFIT_GIT_VERSION "unknown")
endif()

if(NOT EXISTS "${GOOFIT_SOURCE_DIR}/extern/cmake_utils/FindThrust.cmake"
   OR NOT EXISTS "${GOOFIT_SOURCE_DIR}/extern/modern_cmake/ModernCMakeUtils.cmake"
   OR NOT EXISTS "${GOOFIT_SOURCE_DIR}/extern/CLI11/CMakeLists.txt"
   OR NOT EXISTS "${GOOFIT_SOURCE_DIR}/extern/catch2/CMakeLists.txt"
   OR NOT EXISTS "${GOOFIT_SOURCE_DIR}/extern/Eigen/CMakeLists.txt"
   OR NOT EXISTS "${GOOFIT_SOURCE_DIR}/extern/FeatureDetector/CMakeLists.txt"
   OR NOT EXISTS "${GOOFIT_SOURCE_DIR}/extern/MCBooster/CMakeLists.txt"
   OR NOT EXISTS "${GOOFIT_SOURCE_DIR}/extern/thrust/README.md")
  message(
    FATAL_ERROR
      "The submodules were not downloaded! GOOFIT_SUBMODULE was turned off or failed. Please update submodules and try again."
  )
endif()

set(CMAKE_MODULE_PATH "${GOOFIT_SOURCE_DIR}/extern/cmake_utils" ${CMAKE_MODULE_PATH})
include("${GOOFIT_SOURCE_DIR}/extern/modern_cmake/ModernCMakeUtils.cmake")
include("${GOOFIT_SOURCE_DIR}/extern/modern_cmake/ModernCUDA.cmake")

# Add clang-tidy if available
option(GOOFIT_TIDY "Add clang-tidy")
option(GOOFIT_TIDY_FIX "Perform fixes for Clang-Tidy - changes source inplace")
if(GOOFIT_TIDY OR GOOFIT_TIDY_FIX)
  find_program(
    CLANG_TIDY_EXE
    NAMES "clang-tidy"
    DOC "Path to clang-tidy executable")

  if(CLANG_TIDY_EXE)
    if(GOOFIT_TIDY_FIX)
      set(DO_CLANG_TIDY "${CLANG_TIDY_EXE}" "-fix")
    else()
      set(DO_CLANG_TIDY "${CLANG_TIDY_EXE}")
    endif()
  endif()
endif()

# Add CCache if available and if supported
if(NOT GOOFIT_TIDY_FIX)
  find_program(CCACHE_PROGRAM ccache)
  if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    set(CMAKE_CUDA_COMPILER_LAUNCHER "${CCACHE_PROGRAM}") # CMake 3.9+
  endif()
endif()

# Add Sanitizers
set(CMAKE_MODULE_PATH "${GOOFIT_SOURCE_DIR}/extern/sanitizers/cmake" ${CMAKE_MODULE_PATH})
find_package(Sanitizers)

set(GOOFIT_CUDA_OR_GROUPSIZE
    "128"
    CACHE STRING "Overrides the default group distribution for Thrust's transform_reduce")
set(GOOFIT_CUDA_OR_GRAINSIZE
    "7"
    CACHE STRING "Overrides the default grain size for Thrust's transform_reduce")

configure_file("${GOOFIT_SOURCE_DIR}/include/goofit/detail/ThrustOverrideConfig.h.in"
               "${GOOFIT_BINARY_DIR}/include/goofit/detail/ThrustOverrideConfig.h")

set(CMAKE_CXX_STANDARD
    "11"
    CACHE STRING
          "Set a version of C++ standard (11,14,...) to use. Must match ROOT if using ROOT.")
set(CMAKE_CUDA_STANDARD
    "11"
    CACHE STRING "Set a version of C++ standard (11,14,...) to use in CUDA.")

### C++ settings ###
macro(GOOFIT_SETUP_STD)
  set(CMAKE_CXX_EXTENSIONS OFF)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  set(CMAKE_CXX_VISIBILITY_PRESET hidden)
endmacro()

goofit_setup_std()

set(GOOFIT_GCC_WARNINGS
    "-Wall -Wextra -Wno-unknown-pragmas -Wno-long-long -Wno-attributes -Wno-sign-compare -Wno-unused-parameter"
)

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GOOFIT_GCC_WARNINGS}")
  # Helpful but irritating: -Wzero-as-null-pointer-constant -Wsuggest-override
  # no-sign-compare can be removed, but will take some work to clean up
  # Same is true for no-unused-parameter
  # if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.0)
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      "RelWithDebInfo"
      CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
            FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel"
                                               "RelWithDebInfo")
endif()

# Printout time taken to compile
option(GOOFIT_TIME_COMPILE "Print time to compile during compilation" OFF)
if(GOOFIT_TIME_COMPILE)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "time")
endif()

### Options ###
set(DEVICE_LISTING CUDA OMP CPP TBB Auto)
set(HOST_LISTING OMP CPP TBB Auto)
mark_as_advanced(DEVICE_LISTING HOST_LISTING)

set(GOOFIT_DEVICE
    Auto
    CACHE STRING "The compute device, options are ${DEVICE_LISTING}")
set(GOOFIT_HOST
    Auto
    CACHE STRING "The compute device, options are ${HOST_LISTING}")

# Requires CMake 3.3 or greater
if(NOT ${GOOFIT_DEVICE} IN_LIST DEVICE_LISTING)
  message(FATAL_ERROR "You must select a device from ${DEVICE_LISTING}, not ${GOOFIT_DEVICE}")
endif()
if(NOT ${GOOFIT_HOST} IN_LIST HOST_LISTING)
  message(FATAL_ERROR "You must select a host from ${HOST_LISTING}, not ${HOST_DEVICE}")
endif()

# Auto device selection
if(GOOFIT_DEVICE STREQUAL Auto)
  include(CheckLanguage)
  check_language(CUDA)
  if(CUDA_FOUND OR CMAKE_CUDA_COMPILER)
    set(GOOFIT_DEVICE
        CUDA
        CACHE STRING "The compute device, options are ${DEVICE_LISTING}" FORCE)
  else()
    find_package(OpenMP QUIET)
    if(TARGET OpenMP::OpenMP_CXX)
      set(GOOFIT_DEVICE
          OMP
          CACHE STRING "The compute device, options are ${DEVICE_LISTING}" FORCE)
    else()
      set(GOOFIT_DEVICE
          CPP
          CACHE STRING "The compute device, options are ${DEVICE_LISTING}" FORCE)
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

# Nicer filling in for GUI
set_property(CACHE GOOFIT_DEVICE PROPERTY STRINGS ${DEVICE_LISTING})
set_property(CACHE GOOFIT_HOST PROPERTY STRINGS ${HOST_LISTING})

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

if(GOOFIT_DEVICE STREQUAL CUDA)
  set(IS_NOT_CUDA FALSE)
else()
  set(IS_NOT_CUDA TRUE)
endif()

option(GOOFIT_KMATRIX "Activate support for KMatrix (slower CUDA builds)" ${IS_NOT_CUDA})

# This is the common library that holds definitions, etc.
add_library(GooFit_Common INTERFACE)
add_library(GooFit::Common ALIAS GooFit_Common)

target_include_directories(GooFit_Common INTERFACE include)
target_include_directories(GooFit_Common INTERFACE "${GOOFIT_BINARY_DIR}/include")

target_compile_definitions(
  GooFit_Common INTERFACE "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_${GOOFIT_DEVICE}")
target_compile_definitions(GooFit_Common INTERFACE "-DMCBOOSTER_BACKEND=${GOOFIT_DEVICE}")
target_compile_definitions(GooFit_Common
                           INTERFACE "-DTHRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_${GOOFIT_HOST}")
# target_compile_definitions(GooFit_Common INTERFACE _USE_MATH_DEFINES)

#define macro to set stack size (in kb)
if(DEFINED GOOFIT_CUDA_STACKSIZE)
  target_compile_definitions(GooFit_Common INTERFACE CUDA_STACKSIZE=${GOOFIT_CUDA_STACKSIZE})
else()
  target_compile_definitions(GooFit_Common INTERFACE CUDA_STACKSIZE=1024)
endif()

option(GOOFIT_DEBUG "Print debugging messages" OFF)
option(GOOFIT_TRACE "Print messages to trace the behavior of GooFit" OFF)

if(GOOFIT_DEBUG)
  target_compile_definitions(GooFit_Common INTERFACE "-DGOOFIT_DEBUG_FLAG=1")
endif()
if(GOOFIT_TRACE)
  target_compile_definitions(GooFit_Common INTERFACE "-DGOOFIT_TRACE_FLAG=1")
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
configure_file("${GOOFIT_SOURCE_DIR}/include/goofit/detail/Backtrace.h.in"
               "${GOOFIT_BINARY_DIR}/include/goofit/detail/Backtrace.h")

set(GOOFIT_ARCH
    Auto
    CACHE STRING "The GPU Architecture, can be Auto, All, Common, a number, or a name")

option(GOOFIT_MPI "Turn on MPI for goofit" OFF)

if(GOOFIT_MPI)
  find_package(MPI REQUIRED)

  target_link_libraries(GooFit_Common INTERFACE MPI::MPI_CXX)
  target_compile_definitions(GooFit_Common INTERFACE "-DGOOFIT_MPI")

  message(STATUS "MPI found. Use the following to run your program")
  message(
    STATUS
      "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} PROCS ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS"
  )
  message(
    STATUS
      "where PROCS is the number of processors on which to execute the program, EXECUTABLE is the MPI program, and ARGS are the arguments to pass to the MPI program."
  )
endif()

# This macro is designed for external packages to auto-add CUDA
# It is also used to setup up the NEW CUDA support
macro(GOOFIT_OPTIONAL_CUDA)
  if(GOOFIT_DEVICE STREQUAL CUDA)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_EXTENSIONS OFF)

    string(APPEND CMAKE_CUDA_FLAGS " --expt-relaxed-constexpr")

  endif()
endmacro()

if(GOOFIT_DEVICE STREQUAL CUDA)
  # The old FindCUDA method
  goofit_optional_cuda()
  set(CUDA_VERSION ${CMAKE_CUDA_COMPILER_VERSION})

  message(
    STATUS
      "Found CUDA LANGUAGE ${CMAKE_CUDA_COMPILER_VERSION} at ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}"
  )

  set(CUDA_TOOLKIT_ROOT_DIR "${CMAKE_CUDA_COMPILER}")
  get_filename_component(CUDA_TOOLKIT_ROOT_DIR "${CUDA_TOOLKIT_ROOT_DIR}" DIRECTORY)
  get_filename_component(CUDA_TOOLKIT_ROOT_DIR "${CUDA_TOOLKIT_ROOT_DIR}" DIRECTORY)

  include(CUDAUtilities)

  set(cuda_lang "$<COMPILE_LANGUAGE:CUDA>")

  # GooFit isn't popular with nvlink
  # Note that this should be the only nvlink command
  string(APPEND CMAKE_CUDA_FLAGS " -Xnvlink=--disable-warnings")
  string(APPEND CMAKE_CUDA_FLAGS " -Xcompiler=-Wno-attributes")

  set(THRUST_INCLUDE_DIR "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")

  find_library(CUDART_LIB cudart HINTS "${CUDA_TOOLKIT_ROOT_DIR}/lib64"
                                       "${CUDA_TOOLKIT_ROOT_DIR}/lib" "${CUDA_TOOLKIT_ROOT_DIR}")

  add_library(CUDA::Libs IMPORTED INTERFACE)
  set_target_properties(
    CUDA::Libs PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}"
                          INTERFACE_LINK_LIBRARIES "${CUDART_LIB}")
endif()

if(NOT DEFINED cuda_lang)
  set(cuda_lang 0)
endif()

# Config helpers
set(rel_mode "$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>")
set(cxx_lang "$<COMPILE_LANGUAGE:CXX>")
set(cxx_lang_rel "$<AND:${cxx_lang},${rel_mode}>")
set(cuda_lang_rel "$<AND:${cuda_lang},${rel_mode}>")

set(GOOFIT_OPTI
    "-march=native"
    CACHE STRING "compiler flags for optimized builds")
string(REPLACE ";" "," GOOFIT_OPTI_CUDA "${GOOFIT_OPTI}")
target_compile_options(
  GooFit_Common INTERFACE $<${cxx_lang_rel}:$<BUILD_INTERFACE:${GOOFIT_OPTI}>>
                          $<${cuda_lang_rel}:$<BUILD_INTERFACE:-Xcompiler=${GOOFIT_OPTI_CUDA}>>)

find_package(Thrust QUIET)
if(NOT THRUST_FOUND)
  set(THRUST_INCLUDE_DIR "${GOOFIT_SOURCE_DIR}/extern/thrust")
  find_package(Thrust REQUIRED)
endif()

option(GOOFIT_FORCE_LOCAL_THRUST
       "GooFit will use the GitHub version of Thrust by default (CUDA 9 and 10's Thrust is buggy)"
       ON)
if(GOOFIT_FORCE_LOCAL_THRUST)
  message(
    STATUS "Forcing local thrust. GOOFIT_FORCE_LOCAL_THRUST=OFF to disable. Requires CUDA < 10.1")
  message(STATUS "Using CUDA from ${THRUST_INCLUDE_DIRS}")
  target_include_directories(GooFit_Common SYSTEM BEFORE
                             INTERFACE "${GOOFIT_SOURCE_DIR}/extern/thrust")
  if(GOOFIT_DEVICE STREQUAL CUDA AND NOT CMAKE_VERSION VERSION_LESS 3.12)
    target_compile_options(GooFit_Common INTERFACE "SHELL:-isystem ${THRUST_INCLUDE_DIRS}")
  else()
    target_include_directories(GooFit_Common SYSTEM INTERFACE "${THRUST_INCLUDE_DIRS}")
  endif()
else()
  target_include_directories(GooFit_Common SYSTEM INTERFACE "${THRUST_INCLUDE_DIRS}")
  message(STATUS "Using Thrust and CUDA from ${THRUST_INCLUDE_DIRS}")
endif()

if(GOOFIT_DEVICE STREQUAL OMP
   OR GOOFIT_HOST STREQUAL OMP
   OR GOOFIT_DEVICE STREQUAL TBB
   OR GOOFIT_HOST STREQUAL TBB)
  find_package(OpenMP)
  if(NOT TARGET OpenMP::OpenMP_CXX)
    message(FATAL_ERROR "Requested OpenMP or TBB and OpenMP not found")
  endif()

  target_link_libraries(GooFit_Common INTERFACE OpenMP::OpenMP_CXX)
endif()

if(GOOFIT_DEVICE STREQUAL TBB OR GOOFIT_HOST STREQUAL TBB)
  find_package(TBB COMPONENTS tbbmalloc tbbmalloc_proxy tbb_preview)
  target_include_directories(goofit_mt SYSTEM INTERFACE "${TBB_INCLUDE_DIRS}")
  target_link_libraries(GooFit_Common INTERFACE ${TBB_LIBRARIES})
endif()

# Include directories are not picked up by FindCUDA
option(GOOFIT_CERNROOT "Look for the ROOT library (turned off for Python)" ON)

if(GOOFIT_CERNROOT)
  find_package(ROOT 6 COMPONENTS Minuit)
  if(ROOT_FOUND)
    message(STATUS "Found ROOT ${ROOT_VERSION} at ${ROOT_INCLUDE_DIRS}")
    target_compile_definitions(GooFit_Common INTERFACE "-DROOT_FOUND")
    target_link_libraries(GooFit_Common INTERFACE ROOT::Libraries ROOT::Minuit)
  endif()
endif()

# If either ROOT not found or Minuit2 not included
if(NOT TARGET ROOT::Minuit2)
  add_subdirectory("extern/Minuit2")
  set_target_properties(Minuit2 PROPERTIES FOLDER extern)
  set_target_properties(Minuit2Math PROPERTIES FOLDER extern)
  target_compile_definitions(Minuit2Math PUBLIC MATHCORE_STANDALONE)
  add_library(ROOT::Minuit2 ALIAS Minuit2)
  message(STATUS "Using standalone Minuit2")
endif()
target_link_libraries(GooFit_Common INTERFACE ROOT::Minuit2)

# Adding simple libraries
set(CLI11_TESTING OFF)
add_subdirectory("extern/CLI11")
mark_as_advanced(CLI_CXX_STD CLI_EXAMPLES CLI_SINGLE_FILE CLI_SINGLE_FILE_TESTS CLI_TESTING)
target_link_libraries(GooFit_Common INTERFACE CLI11::CLI11)

add_subdirectory("extern/FeatureDetector")
set_target_properties(FeatureDetector PROPERTIES FOLDER extern)
target_link_libraries(GooFit_Common INTERFACE FeatureDetector)

target_compile_options(FeatureDetector PRIVATE $<BUILD_INTERFACE:${GOOFIT_OPTI}>)

## Format
add_subdirectory("extern/fmt")
set_target_properties(fmt PROPERTIES FOLDER extern)
mark_as_advanced(FMT_CPPFORMAT FMT_DOC FMT_INSTALL FMT_PEDANTIC FMT_TEST FMT_USE_CPP11)
target_link_libraries(GooFit_Common INTERFACE fmt)

add_library(rang INTERFACE)
target_include_directories(rang SYSTEM INTERFACE "${GOOFIT_SOURCE_DIR}/extern/rang/include")
target_link_libraries(GooFit_Common INTERFACE rang)

add_library(Eigen INTERFACE)
target_include_directories(Eigen SYSTEM INTERFACE "${GOOFIT_SOURCE_DIR}/extern/Eigen")
target_link_libraries(GooFit_Common INTERFACE Eigen)

# This if(ROUND_FOUND) wrapper in theory should not be needed
# But for some reason FALSE is not "falsy" to the define command
if(ROOT_FOUND)
  set(GOOFIT_ROOT_FOUND ROOT_FOUND)
endif()
set(GOOFIT_ARCH_FLAGS ARCH_FLAGS)

option(GOOFIT_SPLASH "Show a unicode splash or a small plain text one" ON)

# cmake_policy(SET CMP0069 NEW) (will be set by using latest CMake version policy already)
include(CheckIPOSupported)
check_ipo_supported(RESULT SUPPORTS_IPO_LOCAL)
message(STATUS "Compiler supports IPO: ${SUPPORTS_IPO_LOCAL}")

# Make this universally available (even when using submodule)
set(SUPPORTS_IPO
    "${SUPPORTS_IPO_LOCAL}"
    CACHE INTERNAL "")

function(GOOFIT_ADD_LIBRARY GNAME)
  set(options STATIC SHARED MODULE)
  set(oneValueArgs "")
  set(multiValueArgs "")
  cmake_parse_arguments(GOOFIT_LIB "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if(GOOFIT_LIB_SHARED)
    set(GOOFIT_LIB_TYPE SHARED)
  elseif(GOOFIT_LIB_MODULE)
    set(GOOFIT_LIB_TYPE MODULE)
  else()
    set(GOOFIT_LIB_TYPE STATIC)
  endif()

  if(GOOFIT_DEVICE STREQUAL CUDA)
    add_library(${GNAME} ${GOOFIT_LIB_TYPE} ${GOOFIT_LIB_UNPARSED_ARGUMENTS})
    set_target_properties(${GNAME} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  else()
    foreach(N ${GOOFIT_LIB_UNPARSED_ARGUMENTS})
      get_filename_component(NEXT ${N} EXT)
      if(NEXT STREQUAL ".cu")
        set_source_files_properties(${N} PROPERTIES LANGUAGE CXX)
      endif()
    endforeach()
    add_library(${GNAME} ${GOOFIT_LIB_TYPE} ${GOOFIT_LIB_UNPARSED_ARGUMENTS})
    add_sanitizers(${GNAME})
    target_compile_options(${GNAME} PRIVATE -x c++)
  endif()
  target_link_libraries(${GNAME} PUBLIC GooFit::Common)
  if(NOT ${GNAME} STREQUAL PDFs)
    set_target_properties(${GNAME} PROPERTIES FOLDER core)
  endif()

  set_target_properties(${GNAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${SUPPORTS_IPO})

  source_group("Source Files" REGULAR_EXPRESSION ".*\\.c[ucp]p?")
  source_group("Header Files\\detail" REGULAR_EXPRESSION ".*/detail/.*\\.h")
  source_group("Header Files\\utilities" REGULAR_EXPRESSION ".*/utilities/.*\\.h")
  source_group("Header Files\\lineshapes" REGULAR_EXPRESSION ".*/lineshapes/.*\\.h")
  source_group("Header Files\\resonances" REGULAR_EXPRESSION ".*/resonances/.*\\.h")
  source_group("Source Files\\detail" REGULAR_EXPRESSION ".*/detail/.*\\.c[ucp]p?")
  source_group("Source Files\\utilities" REGULAR_EXPRESSION ".*/utilities/.*\\.c[ucp]p?")
  source_group("Source Files\\lineshapes" REGULAR_EXPRESSION ".*/lineshapes/.*\\.c[ucp]p?")
  source_group("Source Files\\resonances" REGULAR_EXPRESSION ".*/resonances/.*\\.c[ucp]p?")

  if(CLANG_TIDY_EXE)
    set_target_properties(${GNAME} PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
  endif()
endfunction()

function(GOOFIT_ADD_EXECUTABLE NAMEEXE)
  if(GOOFIT_DEVICE STREQUAL CUDA)
    add_executable(${NAMEEXE} ${ARGN})
  else()
    foreach(N ${ARGN})
      get_filename_component(NEXT ${N} EXT)
      if(NEXT STREQUAL ".cu")
        set_source_files_properties(${N} PROPERTIES LANGUAGE CXX)
      endif()
    endforeach()
    add_executable(${NAMEEXE} ${ARGN})
    add_sanitizers(${NAMEEXE})
    target_compile_options(${NAMEEXE} PUBLIC -x c++)
  endif()
  target_link_libraries(${NAMEEXE} PUBLIC GooFit::GooFit)
  set_target_properties(${NAMEEXE} PROPERTIES FOLDER projects)
  source_group("Source Files" REGULAR_EXPRESSION ".*\\.c[uc]")

  set_target_properties(${GNAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${SUPPORTS_IPO})

  if(CLANG_TIDY_EXE)
    set_target_properties(${GNAME} PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
  endif()
endfunction()

# This links a file
function(GOOFIT_ADD_LINK)
  if(MSVC) # Not officially supported, but needed to even configure on Windows
    set(COMM copy)
  else()
    set(COMM create_symlink)
  endif()

  foreach(NAMELINK ${ARGN})
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${NAMELINK}")

      execute_process(
        COMMAND ${CMAKE_COMMAND} -E ${COMM} "${CMAKE_CURRENT_SOURCE_DIR}/${NAMELINK}"
                "${CMAKE_CURRENT_BINARY_DIR}/${NAMELINK}"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

    else()
      message(WARNING "${NAMELINK} does not exist. Not making link.")
    endif()
  endforeach()
endfunction()

# This links all files in the current directory to the build directory
function(GOOFIT_ADD_DIRECTORY)
  file(
    GLOB directory_listing
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    *)
  set(skip_files CMakeLists.txt CMakeFiles Makefile makefile .gitignore .git)

  foreach(NAMELINK ${directory_listing})
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${NAMELINK}/CMakeLists.txt")
      # Pass : this should not link subdirectories that are CMake packages
    elseif("${CMAKE_CURRENT_SOURCE_DIR}/${NAMELINK}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
      # Pass : this should not link the current build directory
    elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${NAMELINK}/CMakeCache.txt")
      # Pass : this should not link subdirectories that are CMake build directories
    elseif(NOT "${NAMELINK}" IN_LIST skip_files)
      goofit_add_link("${NAMELINK}")
    endif()
  endforeach()
endfunction()

# This allows goofit to add packages to the package registry and "short circuit" them if needed
macro(GOOFIT_ADD_PACKAGE PACKAGE_NAME)
  option("GOOFIT_PACKAGE_${PACKAGE_NAME}" "Enable building the ${PACKAGE_NAME} package" ON)
  feature_option("GOOFIT_PACKAGE_${PACKAGE_NAME}")

  goofit_add_directory()

  # Since you are not in the correct git directory,
  # make sure git does not work inside the binary folder
  execute_process(COMMAND ${CMAKE_COMMAND} -E touch "${CMAKE_CURRENT_BINARY_DIR}/.git")

  set(OPTS "${ARGN}")
  if(NOT GOOFIT_PACKAGE_${PACKAGE_NAME})
    return()
  elseif(ROOT IN_LIST OPTS AND NOT ROOT_FOUND)
    message(
      WARNING "Requested package ${PACKAGE_NAME} but ROOT was not found! Not building package.")
    return()
  elseif(OLD_CUDA IN_LIST OPTS AND GOOFIT_DEVICE STREQUAL CUDA)
    message(
      WARNING "Requested package ${PACKAGE_NAME} but NEW_CUDA is enabled! Not building package.")
    return()
  endif()
endmacro()

# This is a prepare-all macro for external packages
macro(GOOFIT_EXTERNAL_PACKAGE)
  goofit_setup_std()
  goofit_optional_cuda()
  goofit_add_directory()
endmacro()

target_include_directories(GooFit_Common INTERFACE extern/MCBooster)
target_include_directories(GooFit_Common INTERFACE extern/generics)

add_subdirectory(src)

add_library(goofit_lib INTERFACE)
target_link_libraries(goofit_lib INTERFACE goofit_base PDFs GooFit::Common)
add_library(GooFit::GooFit ALIAS goofit_lib)

option(GOOFIT_TESTS "Build the goofit tests" ${CUR_PROJ})
if(GOOFIT_TESTS)
  enable_testing()
  add_subdirectory(extern/catch2)
  add_subdirectory(tests)
endif()

option(GOOFIT_EXAMPLES "Build the example programs" ${CUR_PROJ})
if(GOOFIT_EXAMPLES)
  add_subdirectory(examples)
endif()

if(EXISTS work)
  add_subdirectory(work)
endif()

option(GOOFIT_PACKAGES "Build any goofit_* packages found" ${CUR_PROJ})
if(GOOFIT_PACKAGES)
  file(
    GLOB list_of_packages
    RELATIVE ${GOOFIT_SOURCE_DIR}
    goofit_*)
  foreach(d ${list_of_packages})
    add_subdirectory(${d})
  endforeach()
endif()

# Output the current GooFit version and compile time options
# Placed after packages to allow packages to change GooFit cached values
configure_file("${GOOFIT_SOURCE_DIR}/include/goofit/Version.h.in"
               "${GOOFIT_BINARY_DIR}/include/goofit/Version.h")
# This changes more often, so is kept separate to speed up rebuilds
configure_file("${GOOFIT_SOURCE_DIR}/include/goofit/VersionGit.h.in"
               "${GOOFIT_BINARY_DIR}/include/goofit/VersionGit.h")

# This allows GOOFIT_PYTHON to be automatic
set(CMAKE_MODULE_PATH "${GOOFIT_SOURCE_DIR}/extern/pybind11/tools" ${CMAKE_MODULE_PATH})
find_package(PythonLibsNew)

if(PYTHONLIBS_FOUND)
  if(EXISTS "${PYTHON_LIBRARIES}")
    option(GOOFIT_PYTHON "Python bindings for goofit" ${CUR_PROJ})
  else()
    message(STATUS "Found Python, but library location ${PYTHON_LIBRARIES} does not exist")
    option(GOOFIT_PYTHON "Python bindings for goofit" OFF)
  endif()
else()
  message(
    STATUS "Did not find the development files for Python, turning off automatic Python bindings")
  option(GOOFIT_PYTHON "Python bindings for goofit" OFF)
endif()

if(GOOFIT_PYTHON)
  if(PYTHONLIBS_FOUND AND EXISTS "${PYTHON_LIBRARIES}")
    message(STATUS "Found Python at ${PYTHON_LIBRARIES}, building bindings")
  else()
    message(
      FATAL_ERROR
        "Did not find Python development libraries, make sure the development package for Python is installed or turn off Python bindings"
    )
  endif()

  unset(PYTHON_VERSION)
  add_subdirectory("extern/pybind11")

  add_subdirectory(python)
endif()

add_feature_info("GOOFIT_DEVICE=${GOOFIT_DEVICE}" ON "The device to compile for")
add_feature_info("GOOFIT_HOST=${GOOFIT_HOST}" ON "The host system to compile for")

feature_option(GOOFIT_KMATRIX)
feature_option(GOOFIT_TRACE)
feature_option(GOOFIT_DEBUG)
feature_option(GOOFIT_MPI)
feature_option(GOOFIT_SPLASH)
feature_option(GOOFIT_TESTS)
feature_option(GOOFIT_EXAMPLES)
feature_option(GOOFIT_PACKAGES)
feature_option(GOOFIT_PYTHON)
feature_option(GOOFIT_TIDY_FIX)

feature_summary(WHAT ENABLED_FEATURES DISABLED_FEATURES PACKAGES_FOUND)
feature_summary(FILENAME ${GOOFIT_BINARY_DIR}/features.log WHAT ALL)

# Packaging support
set(CPACK_PACKAGE_NAME "GooFit")
set(CPACK_PACKAGE_VENDOR "github.com/GooFit/GooFit")
set(CPACK_PACKAGE_CONTACT "https://${CPACK_PACKAGE_VENDOR}")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "GPU/OpenMP fitting package")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_SOURCE_GENERATOR "TGZ")
# CPack collects *everything* except what's listed here.
set(CPACK_SOURCE_IGNORE_FILES
    /.git
    /dist
    /.*build.*
    /\\\\.DS_Store
    /.*\\\\.egg-info
    /var
    /azure-pipelines.yml
    /.ci
    /.travis.yml
    /examples
    /python/examples
    /extern/CLI11/extern
    /extern/catch2/projects
    /extern/fmt/doc
    /extern/pybind11/docs
    /extern/Eigen/doc
    .swp
    /.all-contributorsrc
    /.appveyor.yml
    /.pre-commit.*yaml
    /*.env)
include(CPack)
