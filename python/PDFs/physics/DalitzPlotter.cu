#include <goofit/PDFs/physics/DalitzPlotter.h>
#include <pybind11/pybind11.h>

using namespace GooFit;
namespace py = pybind11;
using namespace pybind11::literals;

void init_DalitzPlotter(py::module &m) {
    py::class_<DalitzPlotter>(m, "DalitzPlotter")
        .def(py::init<GooPdf *, DalitzPlotPdf *>())
        .def("getNumEvents", &DalitzPlotter::getNumEvents)
        .def("getX", &DalitzPlotter::getX, "event"_a)
        .def("getY", &DalitzPlotter::getY, "event"_a)
        .def("getXval", &DalitzPlotter::getXval, "event"_a)
        .def("getYval", &DalitzPlotter::getYval, "event"_a)
        .def("getZval", &DalitzPlotter::getZval, "event"_a)
        .def("getVal", &DalitzPlotter::getVal, "event"_a, "num"_a = 0)
        .def("getDataSet", &DalitzPlotter::getDataSet);
}
