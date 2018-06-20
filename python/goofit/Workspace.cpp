#include <goofit/Python.h>

#include <goofit/DataSet.h>
#include <goofit/PdfBase.h>
#include <goofit/Workspace.h>

using namespace GooFit;

void init_Workspace(py::module &m) {
    py::class_<Workspace>(m, "Workspace")
        .def(py::init<>())
        .def("add_pdf", &Workspace::add_pdf)
        .def("add_var", &Workspace::add_var)
        .def("add_obs", (void (Workspace::*)(Observable)) & Workspace::add_obs)
        .def("add_obs", (void (Workspace::*)(EventNumber)) & Workspace::add_obs)
        .def("add_dataset", &Workspace::add_dataset)
        .def("pdf", &Workspace::pdf, "name"_a)
        .def("var", &Workspace::var, "name"_a)
        .def("obs", &Workspace::obs, "name"_a)
        .def("dataset", &Workspace::dataset, "i"_a)
        .def("__str__", &Workspace::Print)

        ;
}
