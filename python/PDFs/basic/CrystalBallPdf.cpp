#include <goofit/Python.h>

#include <goofit/PDFs/basic/CrystalBallPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/basic/CrystalBallPdf.h>

using namespace GooFit;

void init_CrystalBallPdf(py::module &m) {
    py::class_<CrystalBallPdf, GooPdf>(m, "CrystalBallPdf")
        .def(py::init<std::string, Observable, Variable, Variable, Variable>(),
             "A Gaussian with a power-law tail on one side.",
             "name"_a,
             "x"_a,
             "m"_a,
             "s"_a,
             "a"_a)

        .def(py::init<std::string, Observable, Variable, Variable, Variable, Variable>(),
             "A Gaussian with a power-law tail on one side.",
             "name"_a,
             "x"_a,
             "m"_a,
             "s"_a,
             "a"_a,
             "power"_a)

        .def_static("help", []() { return HelpPrinter(CrystalBallPdf_docs); });
}
