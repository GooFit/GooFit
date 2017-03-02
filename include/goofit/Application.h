#include "goofit/detail/CLI11.hpp"
#include "goofit/Version.h"
#include "thrust/detail/config/device_system.h"
#include "goofit/detail/rang.hpp"

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
    int gpuDev_ = 0;
    bool show_gpus_;
    bool show_threads_;
    int argc_;
    char** argv_;

public:
    /// Make a new Application
    Application(std::string discription,
            int argc, char** argv) : App(discription), argc_(argc), argv_(argv) {

        #ifdef GOOFIT_MPI
        MPI_Init(&argc_, &argv_);
        #endif

        // Fallthrough is useful for most models of GooFit subcommands
        fallthrough();

        #if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
        add_option("--gpu-dev", gpuDev_, "GPU device to use", true)->group("GooFit");
        add_flag("--show-gpus", show_gpus_, "Show the available GPU devices and exit")->group("GooFit");
        #endif
        add_flag("--goofit-details", show_threads_, "Output system and threading details")->group("GooFit");

        add_config("--config", "config.ini", "An ini file with command line options in it");

        // Reset color on exit
        std::atexit([](){std::cout << rang::style::reset;});
        
    }

    /// Shortcut for the lazy
    Application(int argc, char** argv) : Application("", argc, argv) {}

    /// Get the set GPU device
    int get_device() const {return gpuDev_;}

    /// simple run since argc and argv are stored 
    void run() {
        parse(argc_, argv_);
    }

    /// Gets called in parse
    virtual void pre_callback() override {

        set_device();

        if(show_threads_) {
            std::cout << "GOOFIT: Version " << GOOFIT_VERSION_MAJOR << "." << GOOFIT_VERSION_MINOR << "." << GOOFIT_VERSION_PATCH << std::endl;
            #if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
            std::cout << "CUDA: Device " << get_device() << std::endl;
            cudaDeviceProp devProp;
            cudaGetDeviceProperties(&devProp, gpuDev_);
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
        if(show_gpus_) {
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

    int exit(const CLI::Error &e) {
        std::cout << (e.exit_code==0 ? rang::fg::blue : rang::fg::red);
        int rval = CLI::App::exit(e);
        std::cout << rang::fg::reset;
        return rval;
    }

    /// Call if the application might fork, otherwise automatic
    /// For example, if explicitly using omp
    void set_device() const {

        #if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
        if(gpuDev_ != 0) {
            cudaSetDevice(gpuDev_);
        }
        #endif
    }

    /// Cleanup MPI
    ~Application() {
        #ifdef GOOFIT_MPI
        MPI_Finalize();
        #endif
    }
};

}
