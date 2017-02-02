

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
  
set(GOOFIT_GENERAL_FLAGS "-fPIC")
# -Wl,--no-undefined,--no-allow-shlib-undefined")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(GOOFIT_CXX_FLAGS ${GOOFIT_CXX_FLAGS} "${GOOFIT_GENERAL_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0 -march=native -O3")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(GOOFIT_CXX_FLAGS "${GOOFIT_GENERAL_FLAGS} -xHost -O3 -march=native")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(GOOFIT_CXX_FLAGS "${GOOFIT_GENERAL_FLAGS} -O3")
else()
    set(GOOFIT_CXX_FLAGS "${GOOFIT_GENERAL_FLAGS} -O3")
endif()

find_package(CUDA 7.5 REQUIRED)
if(CUDA_FOUND)
else()
    include_directories(${PROJECT_SOURCE_DIR}/include/goofit/fakecuda)
endif()

find_package(OpenMP)
find_package(TBB)

find_package(Thrust REQUIRED)
find_package(ROOT REQUIRED)

if(CUDA_FOUND)
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -std=c++11;
        -I${GOOFIT_INCLUDE_DIR}; --cudart; -O4;
        --expt-relaxed-constexpr; -ftemplate-backtrace-limit=0;
        --generate-line-info; -Xptxas -fmad=true; -Xptxas -dlcm=cg;
        -Xptxas --opt-level=4 )

    set(CUDA_SEPARABLE_COMPILATION OFF)
    set(CUDA_VERBOSE_BUILD ON)
    
    include(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/FindCudaArch.cmake)

    select_nvcc_arch_flags(NVCC_FLAGS_EXTRA)

    list(APPEND CUDA_NVCC_FLAGS ${NVCC_FLAGS_EXTRA})

    set(GOOFIT_CUDA_OPTIONS -Xcompiler ${OpenMP_CXX_FLAGS} -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA  -DTHRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_OMP -lgomp)

    add_definitions("-DMCBOOSTER_BACKEND=CUDA")
endif()


include_directories(${PROJECT_SOURCE_DIR}/include)

