#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/DataSet.h>

namespace py = pybind11;

template<class DataSetBase = DataSet>
class PyDataSet : public DataSetBase {
    using DataSetBase::DataSetBase;
    void addEventVector(std::vector<fptype>& vals, fptype weight = 1) override {
        PYBIND11_OVERLOAD_PURE(void, DataSetBase, addEventVector, vals, weight);
    };
};


void init_DataSet(py::module &m) {
    py::class_<DataSet, PyDataSet<>>(m, "DataSet")
        .def(py::init<Variable*>())
        .def(py::init<Variable*, std::string>())
        .def(py::init<std::vector<Variable*>&>())
        .def(py::init<std::vector<Variable*>&, std::string>())
        .def("addEvent", (void (DataSet::*)()) &DataSet::addEvent)
        .def("addEvent", (void (DataSet::*)(fptype)) &DataSet::addEvent)
        ;
}
