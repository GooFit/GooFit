#pragma once

#include <CLI/Error.hpp>

#include <fmt/format.h>

namespace GooFit {

/// Thrown when a general error is encountered
struct GeneralError : public CLI::Error {
    template<typename... Args>
	GeneralError(std::string name, Args&&... args)
    : Error("GeneralError", fmt::format(name, std::forward<Args>(args)...), 2) {}
    
};
}
