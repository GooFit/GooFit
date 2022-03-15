#include <fmt/format.h>

#include <goofit/Application.h>
#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/detail/Globals.h>
#include <goofit/Version.h>
#include <goofit/VersionGit.h>
#include <goofit/detail/fenv.h>

#include <thrust/detail/config/device_system.h>

#ifndef _MSC_VER
#include <x86/cpu_x86.h>
#endif

#ifdef GOOFIT_MPI
#include <mpi.h>
#endif

#if GOOFIT_ROOT_FOUND
#include <RVersion.h>
#endif

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
#include <cuda_runtime.h>
#endif

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_OMP || THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_TBB
#include <omp.h>
#endif

namespace GooFit {

void signal_handler(int s) {
    std::cout << std::endl << reset << red << bold;
    std::cout << "GooFit: Control-C detected, exiting..." << reset << std::endl;

    cleanup();

    std::exit(1);
}

void print_splash() {
    std::cout << reset << green << "       Welcome to";
#if GOOFIT_SPLASH
    // Just in case, for clang format:
    // clang-format off
    std::string splash = R"raw(
     ██████╗                 ████████╗
    ██╔════╝  █████╗  █████╗ ██╔═════╝  ██╗
    ██║  ███╗██╔══██╗██╔══██╗█████╗██╗██████╗
    ██║   ██║██║  ██║██║  ██║██╔══╝██║╚═██╔═╝
    ╚██████╔╝╚█████╔╝╚█████╔╝██║   ██║  ██║
     ╚═════╝  ╚════╝  ╚════╝ ╚═╝   ╚═╝  ██║
                                   ███████║
                                   ╚══════╝
)raw";
    // clang-format on

    std::cout << reset << dim;
    bool cur_green = false;
    for(int i = 0; i < splash.size(); i++) {
        std::string letter = splash.substr(i, 3);
        bool is_edge
            = letter == "╚" || letter == "╝" || letter == "╗" || letter == "╔" || letter == "║" || letter == "═";
        bool is_block = letter == "█";

        if(is_block && !cur_green) {
            std::cout << reset << green;
            cur_green = true;
        } else if(is_edge && cur_green) {
            std::cout << reset << dim;
            cur_green = false;
        }
        std::cout << splash[i];
        if(splash[i] == '\n')
            std::cout << std::flush;
    }
#else
    std::cout << " GooFit\n";
#endif

    std::cout << reset << std::flush;
}

auto goofit_info_version() -> std::string {
    return fmt::format("GooFit: Version {} ({}) Commit: {}", GOOFIT_VERSION, GOOFIT_TAG, GOOFIT_GIT_VERSION);
}

auto goofit_info_device(int gpuDev_) -> std::string {
    std::string output;
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
    if(gpuDev_ >= 0) {
        output += fmt::format("CUDA {}.{}\n", CUDART_VERSION / 1000, (CUDART_VERSION % 100) / 10.);

        int nDev = 0;
        cudaGetDeviceCount(&nDev);
        output += fmt::format("CUDA: Number of devices: {}\n", nDev);

        if(nDev > 0) {
            cudaDeviceProp devProp;
            cudaGetDeviceProperties(&devProp, gpuDev_);
            output += fmt::format("CUDA: Device {}: {}\n", gpuDev_, devProp.name);

            output += fmt::format("CUDA: Compute {}.{}\n", devProp.major, devProp.minor);
            output += fmt::format("CUDA: Total global memory: {} GB\n", devProp.totalGlobalMem / 1.0e9);
            output += fmt::format("CUDA: Multiprocessors: {}", devProp.multiProcessorCount);

#ifdef GOOFIT_DEBUG_FLAG
            output += fmt::format("\nCUDA: Total amount of shared memory per block: {}\n", devProp.sharedMemPerBlock);
            output += fmt::format("CUDA: Total registers per block: {}\n", devProp.regsPerBlock);
            output += fmt::format("CUDA: Warp size: {}\n", devProp.warpSize);
            output += fmt::format("CUDA: Maximum memory pitch: {}\n", devProp.memPitch);
            output += fmt::format("CUDA: Total amount of constant memory: {}", devProp.totalConstMem);
#endif
        }
    }

#elif THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_OMP
    output += fmt::format("OMP: Number of threads: {}", omp_get_max_threads());
#elif THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_TBB
    output += fmt::format("TBB: Backend selected");
#elif THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CPP
    output += fmt::format("CPP: Single threaded mode");
#endif
    return output;
}

void print_goofit_info(int gpuDev_) {
    std::cout << green;
    std::cout << goofit_info_version() << std::endl;
    std::cout << goofit_info_device(gpuDev_) << std::endl << reset;

#if GOOFIT_ROOT_FOUND
    GOOFIT_INFO("ROOT: {}.{}.{}", ROOT_VERSION_CODE / 65536, ROOT_VERSION_CODE / 256 % 256, ROOT_VERSION_CODE % 256);
#else
    GOOFIT_INFO("ROOT: Not found");
#endif

    // Print out warnings if not fully optimized
    std::cout << red;
    FeatureDetector::cpu_x86::print_warnings();
    std::cout << GooFit::reset << std::flush;

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);

    for(int i = 0; i < deviceCount; i++) {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, i);
        GOOFIT_INFO("CUDA: {} {}: Compute {}.{}; Memory {} GB",
                    i,
                    deviceProp.name,
                    deviceProp.major,
                    deviceProp.minor,
                    deviceProp.totalGlobalMem / pow(1024.0, 3.0));
    }
#endif
}

Application::Application(std::string description, int argc, char **argv)
    : App(description)
    , argc_(argc)
    , argv_(argv) {
#ifdef GOOFIT_MPI
    GOOFIT_TRACE("MPI startup");
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
    if(deviceCount == 0) {
        throw GooFit::Error("Cannot run with no GPUs");
    } else if(deviceCount == 1 && localRank > 1) {
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
    add_option("--gpu-dev", gpuDev_, "GPU device to use")->capture_default_str()->group("GooFit");
#endif
    add_flag_function(
        "--info-only",
        [this](int i) {
            print_splash();
            print_goofit_info(-1);
            throw CLI::Success();
        },
        "Show the available GPU devices and exit")
        ->group("GooFit");
#endif
    auto quiet = add_flag("-q,--quiet", quiet_, "Reduce the verbosity of the Application")->group("GooFit");

    add_flag("--nosplash", splash_, "Do not print a splash")->group("GooFit")->excludes(quiet);

    set_config("--config", "config.ini", "An ini file with command line options in it")->group("GooFit");

    // Reset color on exit (but not control-c)
    std::atexit([]() { std::cout << GooFit::reset; });

#ifndef _MSC_VER
    sigIntHandler.sa_handler = signal_handler;
    sigemptyset(&sigIntHandler.sa_mask);
    sigIntHandler.sa_flags = 0;
    sigaction(SIGINT, &sigIntHandler, nullptr);
#endif
}

void Application::pre_callback() {
    set_device();

    if(!quiet_) {
        if(!splash_)
            print_splash();
        print_goofit_info(gpuDev_);
    }
}

void Application::run() { parse(argc_, argv_); }

void Application::set_device() const {
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
    if(gpuDev_ >= 0) {
        cudaSetDevice(gpuDev_);
    }
#endif
}

auto Application::exit(const CLI::Error &e) -> int {
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

Application::~Application() {
    cleanup();
#ifdef GOOFIT_MPI
    MPI_Finalize();
#endif
}

// This function call is enabled for macOS, too. Will not have an affect for CUDA code.
void Application::set_floating_exceptions() {
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
    GOOFIT_INFO("CUDA does not support floating point exceptions. Please recompile in OMP or CPP mode.");
#else
    feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif
}

auto Application::get_filename(const std::string &input_str, std::string base) const -> std::string {
    // Run from current directory
    if(CLI::ExistingFile(input_str).empty())
        return input_str;

    // Run from a different directory
    std::string prog_name{argv_[0]};
    size_t loc = prog_name.rfind("/");
    if(loc != std::string::npos) {
        std::string cdir = prog_name.substr(0, loc);
        std::string new_input{cdir + "/" + input_str};
        if(CLI::ExistingFile(new_input).empty()) {
            std::cout << "Found file at " << new_input << std::endl;
            return new_input;
        }
    }

    // If all else fails, try to get file from source directory (if base is specified)
    if(base.size() > 0) {
        std::string new_input = std::string(GOOFIT_SOURCE_DIR) + "/" + base + "/" + input_str;
        if(CLI::ExistingFile(new_input).empty()) {
            std::cout << "Found file at " << new_input << std::endl;
            return new_input;
        }
    }

    // If in GooFit as a package, check the goofit parent dir
    std::string new_input = std::string(GOOFIT_SOURCE_DIR) + "/../" + input_str;
    if(CLI::ExistingFile(new_input).empty()) {
        std::cout << "Found file at " << new_input << std::endl;
        return new_input;
    }

    // Could not find file
    throw GooFit::FileError(input_str);
}

} // namespace GooFit
