#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/UnbinnedDataSet.h>

namespace py = pybind11;
using namespace GooFit;

void init_UnbinnedDataSet(py::module &m) {
    py::class_<UnbinnedDataSet, DataSet>(m, "UnbinnedDataSet")
        .def("__init__", [](UnbinnedDataSet &instance, py::args args, py::kwargs kwargs){
            std::string name;
            std::vector<Variable*> vars;
            for(auto arg : args)
                vars.push_back(arg.cast<Variable*>());
            if(kwargs.contains("name"))
                name = kwargs["name"].cast<std::string>();
            new (&instance) UnbinnedDataSet(vars, name);
        });
}
