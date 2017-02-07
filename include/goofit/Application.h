#include "goofit/detail/CLI11.hpp"
#include "goofit/Version.h"
#include "thrust/detail/config/device_system.h"

#ifdef GOOFIT_MPI
#include <mpi.h>
#endif

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
#include <cuda.h>
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
    int argc;
    char** argv;

public:
    /// Make a new Application
    Application(std::string discription,
            int argc, char** argv) : App(discription), argc(argc), argv(argv) {

        #ifdef GOOFIT_MPI
        MPI_Init(&argc, &argv);
        #endif

        #if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
        add_option("--gpu-dev", gpuDev, "GPU device to use", true);
        add_flag("--show-gpus", show_gpus, "Show the available GPU devices and exit");
        #endif
        add_flag("--goofit-details", show_threads, "Output system and threading details");
    }

    /// Shortcut for the lazy
    Application(int argc, char** argv) : Application("", argc, argv) {}

    /// Get the set GPU device
    int get_device() const {return gpuDev;}

    /// simpler run since argc and argv are stored 
    void run() {
        CLI::App::run(argc, argv);
    }

    /// Gets called in parse
    virtual void pre_callback() override {

        set_device();

        if(show_threads) {
            std::cout << "GOOFIT: Version " << GOOFIT_VERSION_MAJOR << "." << GOOFIT_VERSION_MINOR << "." << GOOFIT_VERSION_PATCH << std::endl;
            #if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
            std::cout << "CUDA: Device " << get_device() << std::endl;
            cudaDeviceProp devProp;
            cudaGetDeviceProperties(&devProp, gpuDev);
            std::cout << "CUDA: Compute " << devProp.major << "." << devProp.minor << std::endl;
            std::cout << "CUDA: Total global memory: " <<  devProp.totalGlobalMem / 1.0e9 << "GB" << std::endl;
            std::cout << "CUDA: Multiprocessors: "<< devProp.multiProcessorCount << std::endl;
            std::cout << "CUDA: Total amount of shared memory per block: " << devProp.sharedMemPerBlock << std::endl;
            std::cout << "CUDA: Total registers per block: " << devProp.regsPerBlock << std::endl;
            std::cout << "CUDA: Warp size: " << devProp.warpSize << std::endl;
            std::cout << "CUDA: Maximum memory pitch: " << devProp.memPitch << std::endl;
            std::cout << "CUDA: Total amount of constant memory: " << devProp.totalConstMem << std::endl;

            #elif THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_OMP
            std::cout << "OMP: Number of threads: " << omp_get_max_threads << std::endl;
            #endif
        }

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
            throw GooFit::Success();
        }
        #endif
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

    /// CLeanup MPI
    ~Application() {
        #ifdef GOOFIT_MPI
        MPI_Finalize();
        #endif
    }
};

}
