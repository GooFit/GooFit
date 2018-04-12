#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <algorithm>

#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>

namespace py = pybind11;
using namespace GooFit;
using namespace pybind11::literals;

void init_UnbinnedDataSet(py::module &m) {
    py::class_<UnbinnedDataSet, DataSet>(m, "UnbinnedDataSet")
        .def(py::init([](py::args args, py::kwargs kwargs) {
            std::string name;
            std::vector<Observable> vars;
            for(auto arg : args)
                vars.push_back(arg.cast<Observable>());
            if(kwargs.contains("name"))
                name = kwargs["name"].cast<std::string>();
            return new UnbinnedDataSet(vars, name);
        }))
        .def("__getitem__",
             [&m](const UnbinnedDataSet &instance, py::object value) {
                 auto numpy    = m.import("numpy");
                 auto matrix   = instance.to_matrix<Eigen::Matrix<fptype, Eigen::Dynamic, Eigen::Dynamic>>();
                 auto pymatrix = py::cast(matrix);
                 return pymatrix.attr("__getitem__")(value);
             })
        .def("from_numpy",
             [](UnbinnedDataSet &instance,
                py::array_t<fptype, py::array::c_style | py::array::forcecast> array,
                bool filter) {
                 auto vars = instance.getObservables();
                 for(int j = 0; j < array.shape(1); j++) {
                     for(int i = 0; i < array.shape(0); i++) {
                         vars.at(i).setValue(array.at(i, j));
                     }
                     if(!filter
                        || std::all_of(std::begin(vars), std::end(vars), [](Observable var) { return bool(var); }))
                         instance.addEvent();
                 }
             },
             R"raw(
            Convert a numpy array into data. Optional filter=True argument will remove values out of range.

            This is deprecated in favor of to/from matrix.
            )raw",
             "array"_a,
             "filter"_a = false)
        .def("to_numpy",
             [](UnbinnedDataSet &instance) {
                 size_t cols = instance.getObservables().size();
                 size_t rows = instance.getNumEvents();
                 py::array_t<fptype> result{{cols, rows}};

                 for(int i = 0; i < cols; i++)
                     for(int j = 0; j < rows; j++)
                         result.mutable_at(i, j) = instance.getValue(instance.getObservables().at(i), j);

                 return result;
             })
        .def("to_matrix", &UnbinnedDataSet::to_matrix<Eigen::Matrix<fptype, Eigen::Dynamic, Eigen::Dynamic>>)
        .def("from_matrix",
             &UnbinnedDataSet::from_matrix<Eigen::Matrix<fptype, Eigen::Dynamic, Eigen::Dynamic>>,
             "Append a matrix to a dataset. The final parameter will be a count if the matrix is missing one row."
             "matrix"_a,
             "filter"_a = false)
        .def("loadEvent", &UnbinnedDataSet::loadEvent, "Load an event into the observables", "event_number"_a);
}
