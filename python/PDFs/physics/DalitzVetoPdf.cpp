#include <goofit/Python.h>

#include <goofit/PDFs/physics/Amp3Body_TD.h>
#include <goofit/PDFs/physics/DalitzVetoPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/physics/DalitzVetoPdf.h>

using namespace GooFit;

void init_DalitzVetoPdf(py::module &m) {
    py::class_<DalitzVetoPdf, GooPdf>(m, "DalitzVetoPdf")
        .def_static("help", []() { return HelpPrinter(DalitzVetoPdf_docs); })

        ;
}
