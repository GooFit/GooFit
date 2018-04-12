#include <pybind11/pybind11.h>

#include <goofit/PDFs/physics/DalitzPlotHelpers.h>

using namespace GooFit;
namespace py = pybind11;

template <typename Iterator>

void init_strided_range(py::module &m) {
    py::class_<strided_range>(m, "strided_range").def
}
