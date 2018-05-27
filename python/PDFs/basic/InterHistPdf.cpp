#include <goofit/Python.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/basic/InterHistPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;

void init_InterHistPdf(py::module &m) {
    py::class_<InterHistPdf, GooPdf>(m, "InterHistPdf")
        .def(py::init<std::string, BinnedDataSet *, std::vector<Variable>, std::vector<Observable>>(),
             "name",
             "x",
             "params",
             "obses",
             py::keep_alive<1, 3>());
}
