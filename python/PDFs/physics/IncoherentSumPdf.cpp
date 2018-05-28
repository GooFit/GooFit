#include <goofit/Python.h>

#include <goofit/PDFs/physics/IncoherentSumPdf.h>
#include <goofit/PDFs/physics/TddpPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/physics/IncoherentSumPdf.h>

using namespace GooFit;

void init_IncoherentSumPdf(py::module &m) {
    py::class_<IncoherentSumPdf, GooPdf>(m, "IncoherentSumPdf")
        .def(py::init<std::string, Observable, Observable, EventNumber, DecayInfo3, GooPdf *>(),
             IncoherentSumPdf_docs.c_str(),
             "name"_a,
             "m12"_a,
             "m13"_a,
             "eventNumber"_a,
             "decay"_a,
             "eff"_a,
             py::keep_alive<1, 6>(),
             py::keep_alive<1, 7>())

        .def_static("help", []() { return HelpPrinter(IncoherentSumPdf_docs); })

        ;
}
