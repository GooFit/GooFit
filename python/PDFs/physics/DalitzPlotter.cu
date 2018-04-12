#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>

#include <goofit/PDFs/physics/DalitzPlotter.h>
#include <goofit/Variable.h>

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
        .def("getDataSet", &DalitzPlotter::getDataSet)
        .def("make2D",
             [](const DalitzPlotter &self) {
                 py::array_t<fptype> result{{self.getM12().getNumBins(), self.getM13().getNumBins()}};

                 // Setting this array to 0 is important, since not all values will be filled!
                 for(size_t i = 0; i < self.getM12().getNumBins(); i++)
                     for(size_t j = 0; j < self.getM13().getNumBins(); j++)
                         result.mutable_at(i, j) = 0;

                 for(size_t j = 0; j < self.getNumEvents(); ++j) {
                     size_t currm12                      = self.getX(j);
                     size_t currm13                      = self.getY(j);
                     double val                          = self.getVal(j);
                     result.mutable_at(currm12, currm13) = val;
                 }
                 return result;
             },
             "Make a 2D array for plotting")
        .def("getExtent",
             [](const DalitzPlotter &self) {
                 std::vector<double> extents = {
                     self.getM12().getLowerLimit(),
                     self.getM12().getUpperLimit(),
                     self.getM13().getLowerLimit(),
                     self.getM13().getUpperLimit(),
                 };
                 return extents;
             },
             "Get the extents for plotting");
}
