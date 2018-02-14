#include <goofit/BinnedDataSet.h>
#include <goofit/Variable.h>
#include <pybind11/pybind11.h>

#include "props.h"

using namespace GooFit;

namespace py = pybind11;

void init_BinnedDataSet(py::module &m) {
    py::class_<BinnedDataSet, DataSet>(m, "BinnedDataSet")
        .def(py::init([](py::args args, py::kwargs kwargs) {
            std::string name;
            std::vector<Observable> vars;
            for(auto arg : args)
                vars.push_back(arg.cast<Observable>());
            if(kwargs.contains("name"))
                name = kwargs["name"].cast<std::string>();
            return new BinnedDataSet(vars, name);
        }))

        .def("getBinCenter", (fptype(BinnedDataSet::*)(size_t, size_t) const) & BinnedDataSet::getBinCenter)
        .def_property_readonly("bin_center",
                               (fptype(BinnedDataSet::*)(size_t, size_t) const) & BinnedDataSet::getBinCenter)
        // clang-format off
            ADD_PROP_RO(bin_number, getBinNumber, BinnedDataSet)
            ADD_PROP_RO(bin_volume, getBinVolume, BinnedDataSet)
            ADD_PROP_RO(bin_error, getBinError, BinnedDataSet)
            ADD_PROP_RO(num_bins, getNumBins, BinnedDataSet)
            ADD_PROP_RO(num_events, getNumEvents, BinnedDataSet);
    // clang-format on
}
