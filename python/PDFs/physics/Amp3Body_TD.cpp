#include <goofit/Python.h>

#include <pybind11/stl.h>

#include <goofit/PDFs/physics/Amp3Body_TD.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/MixingTimeResolution.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/physics/Amp3Body_TD.h>

using namespace GooFit;

void init_Amp3Body_TD(py::module &m) {
    py::class_<Amp3Body_TD, Amp3BodyBase> cls(m, "Amp3Body_TD");
    cls.def(py::init<std::string,
                     Observable,
                     Observable,
                     Observable,
                     Observable,
                     EventNumber,
                     DecayInfo3t,
                     MixingTimeResolution *,
                     GooPdf *>(),
            "name"_a,
            "dtime"_a,
            "sigmat"_a,
            "m12"_a,
            "m13"_a,
            "eventNumber"_a,
            "decay"_a,
            "r"_a,
            "eff"_a,
            py::keep_alive<1, 8>(),
            py::keep_alive<1, 9>(),
            py::keep_alive<1, 10>())
        .def(py::init<std::string,
                      Observable,
                      Observable,
                      Observable,
                      Observable,
                      EventNumber,
                      DecayInfo3t,
                      MixingTimeResolution *,
                      GooPdf *,
                      Observable *,
                      Observable *>(),
             "name"_a,
             "dtime"_a,
             "sigmat"_a,
             "m12"_a,
             "m13"_a,
             "eventNumber"_a,
             "decay"_a,
             "r"_a,
             "eff"_a,
             "mistag"_a,
             "charmtag"_a,
             py::keep_alive<1, 7>(),
             py::keep_alive<1, 8>(),
             py::keep_alive<1, 9>(),
             py::keep_alive<1, 10>(),
             py::keep_alive<1, 11>(),
             py::keep_alive<1, 12>())
        .def(py::init<std::string,
                      Observable,
                      Observable,
                      Observable,
                      Observable,
                      EventNumber,
                      DecayInfo3t,
                      std::vector<MixingTimeResolution *> &,
                      GooPdf *,
                      Observable>(),
             "name"_a,
             "dtime"_a,
             "sigmat"_a,
             "m12"_a,
             "m13"_a,
             "eventNumber"_a,
             "decay"_a,
             "r"_a,
             "eff"_a,
             "md0"_a,
             py::keep_alive<1, 7>(),
             py::keep_alive<1, 9>(),
             py::keep_alive<1, 10>())
        .def(py::init<std::string,
                      Observable,
                      Observable,
                      Observable,
                      Observable,
                      EventNumber,
                      DecayInfo3t,
                      std::vector<MixingTimeResolution *> &,
                      GooPdf *,
                      Observable,
                      Observable *,
                      Observable *>(),
             "name"_a,
             "dtime"_a,
             "sigmat"_a,
             "m12"_a,
             "m13"_a,
             "eventNumber"_a,
             "decay"_a,
             "r"_a,
             "eff"_a,
             "md0"_a,
             "mistag"_a,
             "charmtag"_a,
             py::keep_alive<1, 7>(),
             py::keep_alive<1, 8>(),
             py::keep_alive<1, 9>(),
             py::keep_alive<1, 10>(),
             py::keep_alive<1, 11>(),
             py::keep_alive<1, 12>())

        .def("setDataSize", &Amp3Body_TD::setDataSize, "dataSize"_a, "evtSize"_a = 5, "offset"_a = 0)
        .def_static("resetCacheCounter", &Amp3Body_TD::resetCacheCounter)
        .def("getD0Fraction", &Amp3Body_TD::getD0Fraction)
        .def("setD0Fraction", &Amp3Body_TD::setD0Fraction, "d0fraction"_a)
        .def("getFractions",
             &Amp3Body_TD::getFractions,
             "Using the current dataset, return the cached fit fraction values")
        .def_static("help", []() { return HelpPrinter(Amp3Body_TD_docs); })

        ;

    m.attr("TddpPdf") = cls;
}
