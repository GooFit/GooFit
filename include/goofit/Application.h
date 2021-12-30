#pragma once

#include <csignal>

#include <CLI/CLI.hpp>

#include <goofit/Color.h>
#include <goofit/Error.h>
#include <goofit/Log.h>

#define GOOFIT_PARSE(app, ...)                                                                                         \
    try {                                                                                                              \
        app.run();                                                                                                     \
    } catch(const GooFit::ParseError &e) {                                                                             \
        return app.exit(e);                                                                                            \
    }

namespace GooFit {

/// will call the correct exit func, no unwinding of the stack though
void signal_handler(int s);

// Importing into the GooFit namespace the main classes from CLI11
using CLI::App;
using CLI::ExistingDirectory;
using CLI::ExistingFile;
using CLI::ExitCodes;
using CLI::FileError;
using CLI::NonexistentPath;
using CLI::Option;
using CLI::ParseError;
using CLI::Success;

/// Optional print for splash
/// Original: Slant Relief from http://patorjk.com/software/taag/#p=testall&f=Wavy&t=GooFit (tightened a bit)
/// New: Block letters
void print_splash();

/// Print out information about GooFit
void print_goofit_info(int gpuDev_ = 0);

/// Get the version of GooFit as a string
auto goofit_info_version() -> std::string;

/// Get the device information
auto goofit_info_device(int gpuDev_) -> std::string;

class Application : public CLI::App {
  protected:
    int gpuDev_ = 0;
    bool quiet_;
    bool splash_;
    int argc_;
    char **argv_;

    /// Handle control-c codes
    struct sigaction sigIntHandler;

  public:
    /// Make a new Application
    Application(std::string description, int argc, char **argv);

    /// Shortcut for the lazy
    Application(int argc, char **argv)
        : Application("", argc, argv) {}

    /// Get the set GPU device
    auto get_device() const -> int { return gpuDev_; }

    /// simple run since argc and argv are stored
    void run();

    /// Gets called in parse
    void pre_callback() override;

    /// Exit
    auto exit(const CLI::Error &e) -> int;

    /// Call if the application might fork, otherwise automatic
    /// For example, if explicitly using omp
    void set_device() const;

    /// Set floating point errors to throw exceptions instead of NaNs
    /// Does not work with CUDA mode.
    static void set_floating_exceptions();

    /// Cleanup MPI if needed
    ~Application() override;

    /// Get a file from the current directory, looks up one and in the true current directory
    /// Base gives a relative path from the source directory
    auto get_filename(const std::string &input_str, std::string base = "") const -> std::string;
};
} // namespace GooFit
