#include <goofit/PdfBase.h>
#include <goofit/Python.h>

#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/detail/Abort.h>

#include <iostream>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace GooFit;

void init_abort(py::module &m) { m.def("abort", &GooFit::abort, "file"_a, "line"_a, "reason"_a, "Pdf"_a); }
