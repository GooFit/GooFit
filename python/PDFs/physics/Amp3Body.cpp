#include <goofit/Python.h>

#include <pybind11/complex.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/DalitzPlotPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/physics/Amp3Body.h>

using namespace GooFit;

void init_Amp3Body(py::module &m) {
    py::class_<Amp3Body, Amp3BodyBase> cls(m, "Amp3Body");
    cls.def(py::init<std::string, Observable, Observable, EventNumber, DecayInfo3, GooPdf *>(),
            Amp3Body_docs.c_str(),
            "name"_a,
            "m12"_a,
            "m13"_a,
            "eventNumber"_a,
            "decay"_a,
            "eff"_a,
            py::keep_alive<1, 6>(), // Important to keep decay alive, to keep PDFs alive
            py::keep_alive<1, 7>())
        .def("setDataSize", &DalitzPlotPdf::setDataSize, "dataSize"_a, "evtSize"_a = 3)
        .def("getCachedWave", &DalitzPlotPdf::getCachedWave, "i"_a)
        .def("sumCachedWave",
             [](DalitzPlotPdf &self, size_t i) { return std::complex<fptype>(self.sumCachedWave(i)); },
             "i"_a)
        .def("fit_fractions",
             &DalitzPlotPdf::fit_fractions,
             "Using the current dataset, return the cached fit fraction values")
        .def_static("help", []() { return HelpPrinter(Amp3Body_docs); })

        ;

    m.attr("DalitzPlotPdf") = cls;
}
