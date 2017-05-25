#include <pybind11/pybind11.h>
#include <goofit/Variable.h>
#include <goofit/BinnedDataSet.h>

using namespace GooFit;

namespace py = pybind11;

void init_BinnedDataSet(py::module &m) {
    py::class_<BinnedDataSet, DataSet>(m, "BinnedDataSet")
        .def(py::init<Variable*>())
        .def(py::init<Variable*, std::string>())
        .def(py::init<std::vector<Variable*>&>())
        .def(py::init<std::vector<Variable*>&, std::string>())
        .def("getBinCenter", (fptype (BinnedDataSet::*) (size_t, size_t) const) &BinnedDataSet::getBinCenter)
        .def("getBinNumber", &BinnedDataSet::getBinNumber)
        .def("getBinVolume", &BinnedDataSet::getBinVolume)
        .def("getBinError", &BinnedDataSet::getBinError)
        .def("getNumBins", &BinnedDataSet::getNumBins)
        .def("getNumEvents", &BinnedDataSet::getNumEvents)
        ;
}
