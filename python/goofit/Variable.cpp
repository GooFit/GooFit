#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "props.h"

#include <fmt/format.h>
#include <sstream>
#include <string>

#include <goofit/Variable.h>

namespace py = pybind11;
using namespace fmt::literals;
using namespace GooFit;

void init_Variable(py::module &m) {
    py::class_<Indexable>(m, "Indexable")
        .def_property_readonly("name", &Indexable::getName)
        // clang-format off
        ADD_PROP(value, getValue, setValue, Indexable)
        ADD_PROP(index, getIndex, setIndex, Indexable)
        ADD_PROP(upperlimit, getUpperLimit, setUpperLimit, Indexable)
        ADD_PROP(lowerlimit, getLowerLimit, setLowerLimit, Indexable)
        // clang-format on
        .def("__repr__", [](const Indexable &v) { return "<Indexable: {}>"_format(v.getName()); })
        .def("__bool__", &Indexable::operator bool);

    py::class_<Observable, Indexable>(m, "Observable")
        .def(py::init<std::string, fptype, fptype>()) // , "name"_a, "min"_a, "max"_a)
        .def_property_readonly("bin_size", &Observable::getBinSize)
        // clang-format off
        ADD_PROP(numbins, getNumBins, setNumBins, Observable)
        // clang-format on
        .def_property_readonly("bin_centers",
                               [](Observable &self) {
                                   std::vector<double> centers(self.getNumBins());
                                   for(size_t i = 0; i < self.getNumBins(); i++)
                                       centers[i] = self.getLowerLimit() + (i + 0.5) * self.getBinSize();
                                   return centers;
                                   ;
                               })
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
        // clang-format off
        ADD_PROP(error, getError, setError, Variable)
        ADD_PROP(fixed, IsFixed, setFixed, Variable)
        ADD_PROP_RO(fitterIndex, getFitterIndex, Variable)
        ADD_PROP_WO(blind, setBlind, Variable)
        // clang-format on
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
