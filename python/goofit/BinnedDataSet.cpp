#include <pybind11/pybind11.h>
#include <goofit/Variable.h>
#include <goofit/BinnedDataSet.h>

using namespace GooFit;

namespace py = pybind11;

void init_BinnedDataSet(py::module &m) {
    py::class_<BinnedDataSet, DataSet>(m, "BinnedDataSet")
        .def("__init__", [](BinnedDataSet &instance, py::args args, py::kwargs kwargs){
            std::string name;
            std::vector<Variable*> vars;
            for(auto arg : args)
                vars.push_back(arg.cast<Variable*>());
            if(kwargs.contains("name"))
                name = kwargs["name"].cast<std::string>();
            new (&instance) BinnedDataSet(vars, name);
        })
        .def("getBinCenter", (fptype(BinnedDataSet::*)(size_t, size_t) const) & BinnedDataSet::getBinCenter)
        .def("getBinNumber", &BinnedDataSet::getBinNumber)
        .def("getBinVolume", &BinnedDataSet::getBinVolume)
        .def("getBinError", &BinnedDataSet::getBinError)
        .def("getNumBins", &BinnedDataSet::getNumBins)
        .def("getNumEvents", &BinnedDataSet::getNumEvents);
}
