#include <pybind11/pybind11.h>
#include <goofit/Variable.h>

#include <fmt/format.h>
#include <string>
#include <sstream>

namespace py = pybind11;
using namespace fmt::literals;
using namespace GooFit;

void init_Variable(py::module &m) {
    py::class_<Indexable>(m, "Indexable")
        .def(py::init<std::string>())
        .def(py::init<std::string, fptype>())
        .def_property_readonly("name", &Indexable::getName)
        .def_property("value", &Indexable::getValue, &Indexable::setValue)
        .def_property("index", &Indexable::getIndex, &Indexable::setIndex)
        .def_property_readonly("fitterIndex", &Indexable::getFitterIndex);

    py::class_<Variable, Indexable>(m, "Variable")
        .def(py::init<std::string, fptype>())
        .def(py::init<std::string, fptype, fptype>())
        .def(py::init<std::string, fptype, fptype, fptype>())
        .def(py::init<std::string, fptype, fptype, fptype, fptype>())
        .def_property("error", &Variable::getError, &Variable::setError)
        .def_property("upperlimit", &Variable::getUpperLimit, &Variable::setUpperLimit)
        .def_property("lowerlimit", &Variable::getLowerLimit, &Variable::setLowerLimit)
        .def_property("getNumBins", &Variable::getNumBins, &Variable::setNumBins)
        .def_property("fixed", &Variable::IsFixed, &Variable::setFixed)
        //.def_property("blind", &Variable::b, &Variable::setBlind)
        .def("__repr__", [](const Variable &v) { return "<Variable: {}>"_format(v.getName()); })
        .def("__str__", [](const Variable &v) {
            std::stringstream os;
            os << v;
            return os.str();
        })
        .def("__bool__", &Variable::operator bool)
    ;
}
