#pragma once

#include <catch.hpp>

#include <goofit/PDFs/GooPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;

std::string capture_stdout(std::function<void()> &&func);

inline Approx operator"" _a(long double val) { return Approx(val); }

inline Approx operator"" _a(unsigned long long val) { return Approx(val); }
