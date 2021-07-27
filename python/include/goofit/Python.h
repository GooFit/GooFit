#pragma once

#ifdef PYBIND11_VERSION_MAJOR
#error "This must be the first include in your PDF!"
#endif

#include <pybind11/pybind11.h>

#include <string>

namespace py = pybind11;
using namespace pybind11::literals;

struct HelpPrinter {
    std::string help_str;
    auto getHelp() const -> std::string { return help_str; }
    HelpPrinter(std::string input)
        : help_str(input) {}
};
