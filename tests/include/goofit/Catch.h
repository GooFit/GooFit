#pragma once

#include <catch2/catch.hpp>

#include <goofit/PDFs/GooPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
using namespace Catch::literals;

std::string capture_stdout(std::function<void()> &&func);
