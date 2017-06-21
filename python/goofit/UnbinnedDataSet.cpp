#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/UnbinnedDataSet.h>

namespace py = pybind11;
using namespace GooFit;

void init_UnbinnedDataSet(py::module &m) {
    py::class_<UnbinnedDataSet, DataSet>(m, "UnbinnedDataSet")
        .def(py::init<Variable *>())
        .def(py::init<Variable *, std::string>())
        .def(py::init<std::vector<Variable *> &>())
        .def(py::init<std::vector<Variable *> &, std::string>());
}
