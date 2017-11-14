#include <goofit/Variable.h>
#include <pybind11/pybind11.h>

#include <fmt/format.h>
#include <sstream>
#include <string>

namespace py = pybind11;
using namespace fmt::literals;
using namespace GooFit;

void init_Variable(py::module &m) {
    py::class_<Indexable>(m, "Indexable")
        .def_property_readonly("name", &Indexable::getName)
        .def_property("value", &Indexable::getValue, &Indexable::setValue)
        .def_property("index", &Indexable::getIndex, &Indexable::setIndex)
        .def_property("upperlimit", &Indexable::getUpperLimit, &Indexable::setUpperLimit)
        .def_property("lowerlimit", &Indexable::getLowerLimit, &Indexable::setLowerLimit)
        .def("__repr__", [](const Indexable &v) { return "<Indexable: {}>"_format(v.getName()); })
        .def("__bool__", &Indexable::operator bool);

    py::class_<Observable, Indexable>(m, "Observable")
        .def(py::init<std::string, fptype, fptype>()) // , "name"_a, "min"_a, "max"_a)
        .def_property("numbins", &Observable::getNumBins, &Observable::setNumBins)
        .def("__repr__", [](const Observable &v) { return "<Observable: {}>"_format(v.getName()); })
        .def("__str__", [](const Observable &v) {
            std::stringstream os;
            os << v;
            return os.str();
        });

    py::class_<Variable, Indexable>(m, "Variable")
        .def(py::init<std::string, fptype>())
        .def(py::init<std::string, fptype, fptype, fptype>())
        .def(py::init<std::string, fptype, fptype, fptype, fptype>())
        .def(py::init<std::string, fptype, fptype, fptype, fptype, bool>())
                // "name"_a, "value"_a, "error"_a, "min"_a, "max"_a, "fixed"_a=false)
        .def_property("error", &Variable::getError, &Variable::setError)
        .def_property("fixed", &Variable::IsFixed, &Variable::setFixed)
        .def_property_readonly("fitterIndex", &Variable::getFitterIndex)
        .def("setBlind", &Variable::setBlind)
        .def("__repr__", [](const Variable &v) { return "<Variable: {}>"_format(v.getName()); })
        .def("__str__", [](const Variable &v) {
            std::stringstream os;
            os << v;
            return os.str();
        });

    py::class_<EventNumber, Observable>(m, "EventNumber")
        .def(py::init<std::string>())
        .def(py::init<std::string, fptype, fptype>()) //, "Create a counting value for event numbers", "name"_a, "min"_a
                                                      //= 0., "max"_a = static_cast<fptype>(EventNumber::maxint))
        .def_property_readonly_static("maxint", [](py::object) { return EventNumber::maxint; })
        .def("__repr__", [](const EventNumber &v) { return "<EventNumber: {}>"_format(v.getName()); });
}
