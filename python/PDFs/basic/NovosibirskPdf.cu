#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/basic/NovosibirskPdf.h>

using namespace GooFit;
namespace py = pybind11;

void init_NovosibirskPdf(py::module &m) {
    py::class_<NovosibirskPdf, GooPdf>(m, "NovosibirskPdf")
        .def(py::init<std::string, Variable*, Variable*, Variable*, Variable*>())
        ;
}






