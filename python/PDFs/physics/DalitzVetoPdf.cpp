#include <goofit/Python.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/physics/Amp3Body_TD.h>
#include <goofit/PDFs/physics/DalitzVetoPdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/physics/DalitzVetoPdf.h>

using namespace GooFit;

void init_DalitzVetoPdf(py::module &m) {
    py::class_<VetoInfo>(m, "VetoInfo")
        .def(py::init<Variable, Variable, DaughterPair>(), "minimum"_a, "maximum"_a, "cyclic_index"_a);
    //.def_readonly("minimum", &VetoInfo::minimum)
    //.def_readonly("maximum", &VetoInfo::maximum)
    //.def_readonly("cyclic_index", &VetoInfo::cyclic_index);

    py::class_<DalitzVetoPdf, GooPdf>(m, "DalitzVetoPdf")
        .def(py::init<std::string,
                      Observable,
                      Observable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      std::vector<VetoInfo>>(),
             DalitzVetoPdf_docs.c_str(),
             "n"_a,
             "_x"_a,
             "_y"_a,
             "motherM"_a,
             "d1m"_a,
             "d2m"_a,
             "d3m"_a,
             "vetos"_a,
             py::keep_alive<1, 8>())
        .def_static("help", []() { return HelpPrinter(DalitzVetoPdf_docs); })

        ;
}
