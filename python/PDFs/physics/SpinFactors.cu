#include <pybind11/pybind11.h>

#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_SpinFactors(py::module &m) {
    py::class_<SpinFactor, GooPdf>(m, "SpinFactor")
        .def(py::init<std::string, SF_4Body, unsigned int, unsigned int, unsigned int, unsigned int>())

        ;
}
