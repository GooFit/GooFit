#pragma once

#include <CLI/Error.hpp>

namespace GooFit {

/// Thrown when a general error is encountered
struct GeneralError : public CLI::Error {
	GeneralError(std::string name)
	: Error("GeneralError", name, 2) {}
};
};
