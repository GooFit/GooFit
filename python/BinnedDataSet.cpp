#include <pybind11/pybind11.h>
#include <goofit/BinnedDataSet.h>

namespace py = pybind11;

void init_BinnedDataSet(py::module &m) {
    py::class_<BinnedDataSet>(m, "BinnedDataSet")
        .def(py::init<Variable*, std::string>())
        .def("getBinCenter", &BinnedDataSet::getBinCenter)
        .def("getBinNumber", &BinnedDataSet::getBinNumber)
        .def("getBinVolume", &BinnedDataSet::getBinVolume)
        .def("getBinError", &BinnedDataSet::getBinError)
        .def("getNumBins", &BinnedDataSet::getNumBins)
        .def("getNumEvents", &BinnedDataSet::getNumEvents)
        ;
}
