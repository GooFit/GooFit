#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/Variable.h>
#include <goofit/DataSet.h>

namespace py = pybind11;
using namespace GooFit;

template<class DataSetBase = DataSet>
class PyDataSet : public DataSetBase {
    using DataSetBase::DataSetBase;
    void addEvent() override { PYBIND11_OVERLOAD_PURE(void, DataSetBase, addEvent); };
};

void init_DataSet(py::module &m) {
    py::class_<DataSet, PyDataSet<>>(m, "DataSet")
        .def("addEvent", (void (DataSet::*)()) & DataSet::addEvent)
        .def("__len__", &DataSet::getNumEvents)
        .def_property_readonly("name", &DataSet::getName)
        .def_property_readonly("variables", &DataSet::getVariables)
        ;
}
