#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <algorithm>

#include <goofit/Variable.h>
#include <goofit/UnbinnedDataSet.h>

namespace py = pybind11;
using namespace GooFit;
using namespace pybind11::literals;

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
        })
        .def("from_numpy", [](UnbinnedDataSet &instance,
                          py::array_t<fptype, py::array::c_style | py::array::forcecast> array,
                          bool filter){
                auto vars = instance.getVariables();
            for(int j=0; j<array.shape(1); j++) {
                for(int i=0; i<array.shape(0); i++) {
                    vars.at(i)->setValue(array.at(i,j));
                }
                if(!filter || std::all_of(std::begin(vars), std::end(vars), [](Variable* var){return bool(*var);}))
                    instance.addEvent();
            }
        }, 
            R"raw(
            Convert a numpy array into data. Optional filter=True argument will remove values out of range.
            )raw",
            "array"_a, "filter"_a=false)
        .def("to_numpy", [](UnbinnedDataSet &instance) {
            size_t cols = instance.getVariables().size();
            size_t rows = instance.getNumEvents();
            py::array_t<fptype> result{{cols, rows}};
            
            for(int i=0; i<cols; i++)
                for(int j=0; j<rows; j++)
                    result.mutable_at(i,j) = instance.getValue(instance.getVariables().at(i), j);
            
            return result;
        })
    ;
}
