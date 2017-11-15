#include <goofit/Application.h>

#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>

#include <iostream>

namespace py = pybind11;
using namespace GooFit;

void init_Application(py::module &m) {
    m.def("print_splash", &print_splash, "Print a splash screen", py::call_guard<py::scoped_ostream_redirect>());
    m.def("print_goofit_info", [](){
        print_splash();
        GOOFIT_INFO("GooFit: Version {} ({}) Commit: {}", GOOFIT_VERSION, GOOFIT_TAG, GOOFIT_GIT_VERSION);
        
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
        int gpuDev_=0;
        cudaDeviceProp devProp;
        cudaGetDeviceProperties(&devProp, gpuDev_);

        GOOFIT_INFO("CUDA: Device {}: {}", gpuDev_, devProp.name);

        GOOFIT_INFO("CUDA: Compute {}.{}", devProp.major, devProp.minor);
        GOOFIT_INFO("CUDA: Total global memory: {} GB", devProp.totalGlobalMem / 1.0e9);
        GOOFIT_INFO("CUDA: Multiprocessors: {}", devProp.multiProcessorCount);

        GOOFIT_DEBUG("CUDA: Total amount of shared memory per block: {}", devProp.sharedMemPerBlock);
        GOOFIT_DEBUG("CUDA: Total registers per block: {}", devProp.regsPerBlock);
        GOOFIT_DEBUG("CUDA: Warp size: {}", devProp.warpSize);
        GOOFIT_DEBUG("CUDA: Maximum memory pitch: {}", devProp.memPitch);
        GOOFIT_DEBUG("CUDA: Total amount of constant memory: {}", devProp.totalConstMem);

#elif THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_OMP
        GOOFIT_INFO("OMP: Number of threads: {}", omp_get_max_threads());
#elif THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_TBB
        GOOFIT_INFO("TBB: Backend selected");
#elif THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CPP
        GOOFIT_INFO("CPP: Single threaded mode");
#endif

        GOOFIT_INFO("ROOT: {}", GOOFIT_ROOT_FOUND ? "Found" : "Not found");

        // Print out warnings if not fully optimized
        std::cout << red;
        FeatureDetector::cpu_x86::print_warnings();
        std::cout << GooFit::reset << std::flush;
        }, "Print GooFit information", py::call_guard<py::scoped_ostream_redirect>());
}
