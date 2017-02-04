#include "goofit/detail/CLI.hpp"
#include "goofit/Version.h"

#ifdef GOOFIT_MPI
#include <mpi.h>
#endif

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA

#endif

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_OMP
#include <omp.h>
#endif


namespace GooFit {

using namespace CLI;

class Application : public CLI::App {
protected:
    int gpuDev = 0;
    bool show_gpus;
    bool show_threads;

public:
    Application(std::string discription) : App(discription) {

        #if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
        add_option("--gpu-dev", gpuDev, "GPU device to use", GooFit::Default);
        add_flag("--show-gpus", show_gpus, "Show the availble GPU devices");
        #endif
        add_flag("--goofit-details", show_threads, "Output system and threading details");
    }

    int get_gpu() const {return gpuDev;}

    /// Overriding parse
    void parse(int argc, char** argv) {
        #ifdef GOOFIT_MPI
        MPI_Init(&argc, &argv);
        #endif

        CLI::App::parse(argc, argv);
    }

    /// Gets called in parse
    virtual void pre_callback() override {
        #if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
        if(show_gpus) {
            int deviceCount;
            cudaGetDeviceCount(&deviceCount);
            for(int i = 0; i < deviceCount; i++) {
                cudaDeviceProp deviceProp;
                cudaGetDeviceProperties(&deviceProp, i);
                std::cout << "CUDA: Device " << i << " has compute capability "
                          << deviceProp.major << "." << deviceProp.minor
                          << std::endl;
            }
        }

        set_device();
        #endif
        if(show_threads) {
            std::cout << "GOOFIT: Version " << GOOFIT_VERSION_MAJOR << "." GOOFIT_VERSION_MINOR << "." GOOFIT_VERSION_PATCH << std::endl;
            #if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
            cudaDeviceProp devProp;
            cudaGetDeviceProperties(&devProp, gpuDev);
            std::cout << "CUDA: Compute " << devProp.major << "." << devProp.minor << std::endl;
            std::cout << "CUDA: Total global memory: " <<  devProp.totalGlobalMem / 1.0e9 << "GB" << std::endl;
            std::cout << "CUDA: Multiprocessors: "<< devProp.multiProcessorCount << std::endl;
            std::cout << "CUDA: Total amount of shared memory per block: " << devProp.sharedMemPerBlock << std::endl;
            std::cout << "CUDA: Total registers per block: ", devProp.regsPerBlock << std::endl;
            std::cout << "CUDA: Warp size: " << devProp.warpSize << std::endl;
            std::cout << "CUDA: Maximum memory pitch: " << devProp.memPitch << std::endl;
            std::cout << "CUDA: Total amount of constant memory: " << devProp.totalConstMem << std::endl;

            #elif THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_OMP
            std::cout << "OMP: Number of threads: " << omp_get_max_threads << std::endl;
            #endif
        }

    }

    /// Call if the application might fork, otherwise automatic
    /// For example, if explicitly using omp
    void set_device() const {

        #if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
        if(gpuDev != 0) {
            cudaSetDevice(gpuDev);
        }
        #endif
    }

    ~Application() {
        #ifdef GOOFIT_MPI
        MPI_Finalize();
        #endif

    }
};

}
