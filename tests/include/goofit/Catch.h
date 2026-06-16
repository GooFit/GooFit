#pragma once

#include <catch2/catch_all.hpp>

#include <goofit/PDFs/GooPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
using namespace Catch::literals;

// Catch2 v3 moved Approx into the Catch namespace; keep the unqualified spelling.
using Catch::Approx;

std::string capture_stdout(std::function<void()> &&func);
