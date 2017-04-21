#include <pybind11/pybind11.h>
#include <goofit/Variable.h>
#include <string>
#include <sstream>

namespace py = pybind11;

void init_Variable(py::module &m) {
    py::class_<Indexable>(m, "Indexable")
        .def(py::init<std::string>())
        .def(py::init<std::string, fptype>())
        .def_readwrite("name", &Variable::name)
        .def_readwrite("value", &Variable::value)
        .def_readwrite("index", &Variable::index)
        ;


    py::class_<Variable, Indexable>(m, "Variable")
        .def(py::init<std::string>())
        .def(py::init<std::string, fptype>())
        .def(py::init<std::string, fptype, fptype>())
        .def(py::init<std::string, fptype, fptype, fptype>())
        .def(py::init<std::string, fptype, fptype, fptype, fptype>())
        .def_readwrite("error", &Variable::error)
        .def_readwrite("upperlimit", &Variable::upperlimit)
        .def_readwrite("numbins", &Variable::numbins)
        .def_readwrite("fixed", &Variable::fixed)
        .def_readwrite("blind", &Variable::blind)
        .def("__repr__", [](const Variable &v){
            return "<Variable: " + v.name + ">";
        })
        .def("__str__", [](const Variable &v){
            std::stringstream os;
            os << v;
            return os.str();
        })
        ;
}
