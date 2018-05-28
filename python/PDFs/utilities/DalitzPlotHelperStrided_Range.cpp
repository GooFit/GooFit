#include <goofit/Python.h>

#include <goofit/PDFs/physics/DalitzPlotHelpers.h>

using namespace GooFit;

template <typename Iterator>

void init_strided_range(py::module &m) {
    py::class_<strided_range>(m, "strided_range");
}
