#pragma once

#include <fmt/format.h>
#include "goofit/Color.h"

#include <iostream>

#define GOOFIT_INFO(...)  std::cout << GooFit::reset << GooFit::blue; fmt::print(__VA_ARGS__); std::cout << GooFit::reset << std::endl
#define GOOFIT_INFO_C(color, ...) std::cout << GooFit::reset << color; fmt::print(__VA_ARGS__); std::cout << GooFit::reset << std::endl
#define GOOFIT_STATUS(...)  std::cout << GooFit::reset << GooFit::magenta; fmt::print(__VA_ARGS__); std::cout << GooFit::reset << std::endl

#ifdef GOOFIT_DEBUG_FLAG
#define GOOFIT_DEBUG(...) std::cout << GooFit::reset << GooFit::cyan << GooFit::bold << "DEBUG: "; fmt::print(__VA_ARGS__); std::cout << GooFit::reset << std::endl
#else
#define GOOFIT_DEBUG(...)
#endif
#ifdef GOOFIT_TRACE_FLAG
#define GOOFIT_TRACE(...) std::cout << GooFit::reset << GooFit::cyan << "TRACE: "; fmt::print(__VA_ARGS__); std::cout << GooFit::reset << std::endl
#else
#define GOOFIT_TRACE(...)
#endif

