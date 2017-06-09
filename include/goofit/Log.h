#pragma once

#include "goofit/Color.h"
#include <fmt/format.h>

#include <iostream>

namespace GooFit {

#define GOOFIT_INFO(...)                                                                                               \
    {                                                                                                                  \
        std::cout << GooFit::reset << GooFit::blue;                                                                    \
        fmt::print(__VA_ARGS__);                                                                                       \
        std::cout << GooFit::reset << std::endl;                                                                       \
    }
#define GOOFIT_INFO_F(...)                                                                                             \
    {                                                                                                                  \
        std::cout << GooFit::reset << GooFit::blue;                                                                    \
        fmt::printf(__VA_ARGS__);                                                                                      \
        std::cout << GooFit::reset << std::endl;                                                                       \
    }
#define GOOFIT_INFO_C(color, ...)                                                                                      \
    {                                                                                                                  \
        std::cout << GooFit::reset << color;                                                                           \
        fmt::print(__VA_ARGS__);                                                                                       \
        std::cout << GooFit::reset << std::endl;                                                                       \
    }
#define GOOFIT_INFO_FC(color, ...)                                                                                     \
    {                                                                                                                  \
        std::cout << GooFit::reset << color;                                                                           \
        fmt::printf(__VA_ARGS__);                                                                                      \
        std::cout << GooFit::reset << std::endl;                                                                       \
    }
#define GOOFIT_STATUS(...)                                                                                             \
    {                                                                                                                  \
        std::cout << GooFit::reset << GooFit::magenta;                                                                 \
        fmt::print(__VA_ARGS__);                                                                                       \
        std::cout << GooFit::reset << std::endl;                                                                       \
    }
#define GOOFIT_STATUS_F(...)                                                                                           \
    {                                                                                                                  \
        std::cout << GooFit::reset << GooFit::magenta;                                                                 \
        fmt::printf(__VA_ARGS__);                                                                                      \
        std::cout << GooFit::reset << std::endl;                                                                       \
    }
#define GOOFIT_WARN(...)                                                                                               \
    {                                                                                                                  \
        std::cout << GooFit::reset << GooFit::yellow << GooFit::bold;                                                  \
        fmt::print(__VA_ARGS__);                                                                                       \
        std::cout << GooFit::reset << std::endl;                                                                       \
    }
#define GOOFIT_WARN_F(...)                                                                                             \
    {                                                                                                                  \
        std::cout << GooFit::reset << GooFit::yellow << GooFit::bold;                                                  \
        fmt::printf(__VA_ARGS__);                                                                                      \
        std::cout << GooFit::reset << std::endl;                                                                       \
    }
#define GOOFIT_ERROR(...)                                                                                              \
    {                                                                                                                  \
        std::cout << GooFit::reset << GooFit::red << GooFit::bold;                                                     \
        fmt::print(__VA_ARGS__);                                                                                       \
        std::cout << GooFit::reset << std::endl;                                                                       \
    }
#define GOOFIT_ERROR_F(...)                                                                                            \
    {                                                                                                                  \
        std::cout << GooFit::reset << GooFit::red << GooFit::bold;                                                     \
        fmt::printf(__VA_ARGS__);                                                                                      \
        std::cout << GooFit::reset << std::endl;                                                                       \
    }

#ifdef GOOFIT_DEBUG_FLAG
#ifndef __CUDA_ARCH__
#define GOOFIT_DEBUG(...)                                                                                              \
    {                                                                                                                  \
        std::cout << GooFit::reset << GooFit::cyan << GooFit::bold << "DEBUG: ";                                       \
        fmt::print(__VA_ARGS__);                                                                                       \
        std::cout << GooFit::reset << std::endl;                                                                       \
    }
#define GOOFIT_DEBUG_F(...)                                                                                            \
    {                                                                                                                  \
        std::cout << GooFit::reset << GooFit::cyan << GooFit::bold << "DEBUG: ";                                       \
        fmt::printf(__VA_ARGS__);                                                                                      \
        std::cout << GooFit::reset << std::endl;                                                                       \
    }
#else
#define GOOFIT_DEBUG(...)
#define GOOFIT_DEBUG_F(...)                                                                                            \
    {                                                                                                                  \
        if(blockId.x == 0 && blockId.y == 0 && threadId.x == 0 && threadId.y == 0) {                                   \
            printf(__VA_ARGS__);                                                                                       \
        }                                                                                                              \
    }
#endif
#else
#define GOOFIT_DEBUG(...)
#define GOOFIT_DEBUG_F(...)
#endif
#ifdef GOOFIT_TRACE_FLAG
#ifndef __CUDA_ARCH__
#define GOOFIT_TRACE(...)                                                                                              \
    {                                                                                                                  \
        std::cout << GooFit::reset << GooFit::cyan << "TRACE: ";                                                       \
        fmt::print(__VA_ARGS__);                                                                                       \
        std::cout << GooFit::reset << std::endl;                                                                       \
    }
#define GOOFIT_TRACE_F(...)                                                                                            \
    {                                                                                                                  \
        std::cout << GooFit::reset << GooFit::cyan << "TRACE: ";                                                       \
        fmt::printf(__VA_ARGS__);                                                                                      \
        std::cout << GooFit::reset << std::endl;                                                                       \
    }
#else
#define GOOFIT_TRACE(...)
#define GOOFIT_TRACE_F(...)                                                                                            \
    {                                                                                                                  \
        if(blockId.x == 0 && blockId.y == 0 && threadId.x == 0 && threadId.y == 0) {                                   \
            printf(__VA_ARGS__);                                                                                       \
        }                                                                                                              \
    }
#endif
#else
#define GOOFIT_TRACE(...)
#define GOOFIT_TRACE_F(...)
#endif

} // namespace GooFit
