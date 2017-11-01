#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/MixingTimeResolution_Aux.h>
#include <goofit/PDFs/physics/TddpPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_TddpPdf(py::module &m) {
    py::class_<TddpPdf, GooPdf>(m, "TddpPdf")
        .def(py::init<std::string,
                      Observable *,
                      Observable *,
                      Observable *,
                      Observable *,
                      EventNumber *,
                      DecayInfo *,
                      MixingTimeResolution *,
                      GooPdf *>(),
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>(),
             py::keep_alive<1, 6>(),
             py::keep_alive<1, 7>(),
             py::keep_alive<1, 8>(),
             py::keep_alive<1, 9>(),
             py::keep_alive<1, 10>())
        .def(py::init<std::string,
                      Observable *,
                      Observable *,
                      Observable *,
                      Observable *,
                      EventNumber *,
                      DecayInfo *,
                      MixingTimeResolution *,
                      GooPdf *,
                      Observable *>(),
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>(),
             py::keep_alive<1, 6>(),
             py::keep_alive<1, 7>(),
             py::keep_alive<1, 8>(),
             py::keep_alive<1, 9>(),
             py::keep_alive<1, 10>(),
             py::keep_alive<1, 11>())
        .def(py::init<std::string,
                      Observable *,
                      Observable *,
                      Observable *,
                      Observable *,
                      EventNumber *,
                      DecayInfo *,
                      std::vector<MixingTimeResolution *> &,
                      GooPdf *,
                      Observable *>(),
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>(),
             py::keep_alive<1, 6>(),
             py::keep_alive<1, 7>(),
             py::keep_alive<1, 8>(),
             py::keep_alive<1, 9>(),
             py::keep_alive<1, 10>(),
             py::keep_alive<1, 11>())
        .def(py::init<std::string,
                      Observable *,
                      Observable *,
                      Observable *,
                      Observable *,
                      EventNumber *,
                      DecayInfo *,
                      std::vector<MixingTimeResolution *> &,
                      GooPdf *,
                      Observable *,
                      Observable *>(),
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>(),
             py::keep_alive<1, 5>(),
             py::keep_alive<1, 6>(),
             py::keep_alive<1, 7>(),
             py::keep_alive<1, 8>(),
             py::keep_alive<1, 9>(),
             py::keep_alive<1, 10>(),
             py::keep_alive<1, 11>(),
             py::keep_alive<1, 12>());
    
    // TODO: Please use annotations and defaults here!
}
