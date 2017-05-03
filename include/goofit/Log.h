#pragma once

#include <fmt/format.h>
#include "goofit/Color.h"

#include <iostream>

#define GOOFIT_INFO(...)  std::cout << GooFit::reset << GooFit::blue; fmt::print(__VA_ARGS__); std::cout << GooFit::reset << std::endl
#define GOOFIT_INFO_C(color, ...) std::cout << GooFit::reset << color; fmt::print(__VA_ARGS__); std::cout << GooFit::reset << std::endl
#define GOOFIT_STATUS(...)  std::cout << GooFit::reset << GooFit::magenta; fmt::print(__VA_ARGS__); std::cout << GooFit::reset << std::endl

#ifdef GOOFIT_DEBUG_FLAG
#ifndef __CUDA_ARCH__
#define GOOFIT_DEBUG(...) std::cout << GooFit::reset << GooFit::cyan << GooFit::bold << "DEBUG: "; fmt::print(__VA_ARGS__); std::cout << GooFit::reset << std::endl
#define GOOFIT_DEBUGF(...) std::cout << GooFit::reset << GooFit::cyan << GooFit::bold << "DEBUG: "; fmt::printf(__VA_ARGS__); std::cout << GooFit::reset << std::endl
#else
#define GOOFIT_DEBUG(...)
#define GOOFIT_DEBUGF(...) if (blockId.x == 0 && blockId.y == 0 && threadId.x == 0 && threadId.y == 0) {printf(__VA_ARGS__);}
#endif
#else
#define GOOFIT_DEBUG(...)
#define GOOFIT_DEBUGF(...)
#endif
#ifdef GOOFIT_TRACE_FLAG
#ifndef __CUDA_ARCH__
#define GOOFIT_TRACE(...) std::cout << GooFit::reset << GooFit::cyan << "TRACE: "; fmt::print(__VA_ARGS__); std::cout << GooFit::reset << std::endl
#define GOOFIT_TRACEF(...) std::cout << GooFit::reset << GooFit::cyan << "TRACE: "; fmt::printf(__VA_ARGS__); std::cout << GooFit::reset << std::endl
#else
#define GOOFIT_TRACE(...)
#define GOOFIT_TRACEF(...) if (blockId.x == 0 && blockId.y == 0 && threadId.x == 0 && threadId.y == 0) {printf(__VA_ARGS__);}
#endif
#else
#define GOOFIT_TRACE(...)
#define GOOFIT_TRACEF(...)
#endif

