#pragma once

#include <CLI/CLI.hpp>
#include "goofit/Version.h"
#include "goofit/Error.h"
#include "goofit/Color.h"
#include "goofit/Log.h"

#include <thrust/detail/config/device_system.h>
#include <x86/cpu_x86.h>

#ifdef GOOFIT_MPI
#include <mpi.h>
#endif

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
#include <cuda.h>
#endif

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_OMP || THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_TBB
#include <omp.h>
#endif

#include <csignal>

namespace GooFit {

void signal_handler(int s) {
    std::cout << std::endl << reset << red << bold;
    std::cout << "GooFit: Control-C detected, exiting..." << reset << std::endl;
    std::exit(1); // will call the correct exit func, no unwinding of the stack though
}

// Importing into the GooFit namespace the main classes from CLI11
using CLI::ParseError;
using CLI::FileError;
using CLI::Success;
using CLI::App;
using CLI::Option;
using CLI::ExistingFile;
using CLI::ExistingDirectory;
using CLI::NonexistentPath;
using CLI::ExitCodes;

class Application : public CLI::App {
  protected:
    int gpuDev_ = 0;
    bool show_gpus_;
    bool quiet_;
    int argc_;
    char **argv_;

    /// Handle control-c codes
    struct sigaction sigIntHandler;

  public:
    /// Make a new Application
    Application(std::string discription, int argc, char **argv)
        : App(discription)
        , argc_(argc)
        , argv_(argv) {
#ifdef GOOFIT_MPI
        MPI_Init(&argc_, &argv_);

        int myId, numProcs;
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myId);

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
        int deviceCount;
        cudaGetDeviceCount(&deviceCount);

        // Get Batch system number of nodes if possible
        auto PBS_NUM_NODES = getenv("PBS_NUM_NODES");

        int nodes = PBS_NUM_NODES == nullptr ? 1 : atoi(PBS_NUM_NODES);

        if(nodes == 0)
            nodes = 1;

        int procsPerNode = numProcs / nodes;
        int localRank    = myId % procsPerNode;

        // Note, this will (probably) be overwritten by gpu-set-device calls...
        if(deviceCount == 1 && localRank > 1) {
            // Multi-process to one GPU!
            gpuDev_ = 0;
        } else if(procsPerNode > 1 && deviceCount > 1) {
            if(localRank <= deviceCount) {
                // setting multiple processes to multiple GPU's
                gpuDev_ = localRank;
            } else {
                // More multiple processes than GPU's, distribute sort of evenly!
                gpuDev_ = localRank % deviceCount;
            }
        } else {
            // multiple GPU's, using one process
            gpuDev_ = 0;
        }

        std::cout << "MPI using CUDA device: " << gpuDev_ << std::endl;
        cudaSetDevice(gpuDev_);
#endif
#endif

        // Fallthrough is useful for most models of GooFit subcommands
        fallthrough();

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
#ifndef GOOFIT_MPI
        add_option("--gpu-dev", gpuDev_, "GPU device to use", true)->group("GooFit");
#endif
        add_flag("--show-gpus", show_gpus_, "Show the available GPU devices and exit")->group("GooFit");
#endif
        add_flag("-q,--quiet", quiet_, "Reduce the verbosity of the Application")->group("GooFit");

        add_config("--config", "config.ini", "An ini file with command line options in it")->group("GooFit");

        // Reset color on exit (but not control-c)
        std::atexit([]() { std::cout << GooFit::reset; });

        sigIntHandler.sa_handler = signal_handler;
        sigemptyset(&sigIntHandler.sa_mask);
        sigIntHandler.sa_flags = 0;
        sigaction(SIGINT, &sigIntHandler, nullptr);
    }

    /// Shortcut for the lazy
    Application(int argc, char **argv)
        : Application("", argc, argv) {}

    /// Get the set GPU device
    int get_device() const { return gpuDev_; }

    /// simple run since argc and argv are stored
    void run() { parse(argc_, argv_); }

    /// Gets called in parse
    void pre_callback() override {
        set_device();

        if(!quiet_) {
            GOOFIT_INFO("GooFit: Version {} ({}) Commit: {}", GOOFIT_VERSION, GOOFIT_TAG, GOOFIT_GIT_VERSION);

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
            cudaDeviceProp devProp;
            cudaGetDeviceProperties(&devProp, gpuDev_);

            GOOFIT_INFO("CUDA: Device {}: {}", get_device(), devProp.name);

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

            // Print out warnings if not fully optimized
            std::cout << red;
            FeatureDetector::cpu_x86::print_warnings();
            std::cout << GooFit::reset << std::flush;
        }

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA

        if(show_gpus_) {
            int deviceCount;
            cudaGetDeviceCount(&deviceCount);

            for(int i = 0; i < deviceCount; i++) {
                cudaDeviceProp deviceProp;
                cudaGetDeviceProperties(&deviceProp, i);
                std::cout << "CUDA: Device " << i << " has compute capability " << deviceProp.major << "."
                          << deviceProp.minor << std::endl;
            }

            throw GooFit::Success();
        }

#endif
    }

    int exit(const CLI::Error &e) {
#ifdef GOOFIT_MPI
        int myId;
        MPI_Comm_rank(MPI_COMM_WORLD, &myId);

        if(myId > 0)
            return e.get_exit_code();

#endif
        std::cout << (e.get_exit_code() == 0 ? blue : red);
        int rval = CLI::App::exit(e);
        std::cout << GooFit::reset;
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

    /// Get a file from the current directory, looks up one and in the true current directory
    /// Base gives a relative path from the source directory
    std::string get_filename(const std::string &input_str, std::string base = "") const {
        // Run from current directory
        if(CLI::ExistingFile(input_str))
            return input_str;

        // Run from a different directory
        std::string prog_name{argv_[0]};
        size_t loc = prog_name.rfind("/");
        if(loc != std::string::npos) {
            std::string cdir = prog_name.substr(0, loc);
            std::string new_input{cdir + "/" + input_str};
            if(CLI::ExistingFile(new_input)) {
                std::cout << "Found file at " << new_input << std::endl;
                return new_input;
            }
        }

        // If all else fails, try to get file from source directory (if base is specified)
        if(base.size() > 0) {
            std::string new_input = std::string(GOOFIT_SOURCE_DIR) + "/" + base + "/" + input_str;
            if(CLI::ExistingFile(new_input)) {
                std::cout << "Found file at " << new_input << std::endl;
                return new_input;
            }
        }

        // Could not find file
        throw GooFit::FileError(input_str);
    }
};
}
