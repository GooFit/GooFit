#include <goofit/Python.h>

#include <goofit/PDFs/physics/Amp3Body_IS.h>
#include <goofit/PDFs/physics/detail/SpecialIncoherentIntegrator.h>
#include <goofit/PDFs/physics/detail/SpecialIncoherentResonanceCalculator.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/physics/Amp3Body_IS.h>

using namespace GooFit;

void init_Amp3Body_IS(py::module &m) {
    py::class_<Amp3Body_IS, GooPdf> cls(m, "Amp3Body_IS");
    cls.def(py::init<std::string, Observable, Observable, EventNumber, DecayInfo3, GooPdf *>(),
            Amp3Body_IS_docs.c_str(),
            "name"_a,
            "m12"_a,
            "m13"_a,
            "eventNumber"_a,
            "decay"_a,
            "eff"_a,
            py::keep_alive<1, 6>(),
            py::keep_alive<1, 7>())

        .def_static("help", []() { return HelpPrinter(Amp3Body_IS_docs); })

        ;

    m.attr("IncoherentSumPdf") = cls;
}
