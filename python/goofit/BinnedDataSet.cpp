#include <goofit/BinnedDataSet.h>
#include <goofit/Variable.h>
#include <pybind11/pybind11.h>

using namespace GooFit;

namespace py = pybind11;

void init_BinnedDataSet(py::module &m) {
    py::class_<BinnedDataSet, DataSet>(m, "BinnedDataSet")
        .def(py::init([](py::args args, py::kwargs kwargs) {
            std::string name;
            std::vector<Variable *> vars;
            for(auto arg : args)
                vars.push_back(arg.cast<Variable *>());
            if(kwargs.contains("name"))
                name = kwargs["name"].cast<std::string>();
            return new BinnedDataSet(vars, name);
        }))
        .def("getBinCenter", (fptype (BinnedDataSet::*)(size_t, size_t) const) & BinnedDataSet::getBinCenter)
        .def("getBinNumber", &BinnedDataSet::getBinNumber)
        .def("getBinVolume", &BinnedDataSet::getBinVolume)
        .def("getBinError", &BinnedDataSet::getBinError)
        .def("getNumBins", &BinnedDataSet::getNumBins)
        .def("getNumEvents", &BinnedDataSet::getNumEvents);
}
