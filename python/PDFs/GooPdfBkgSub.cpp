#include <goofit/Python.h>

#include <pybind11/stl.h>

#include <goofit/PDFs/GooPdfBkgSub.h>
#include <goofit/Variable.h>

using namespace GooFit;

void init_GooPdfBkgSub(py::module &m) {
    py::class_<GooPdfBkgSub, GooPdf>(m, "GooPdfBkgSub")
        .def(py::init<GooPdf, UnbinnedDataSet*, Variable>(),
             "", // Not sure how to deal with generating this. Will ignore for now.
             "basePdf"_a,
             "bkgEvents"_a,
             "bkgFraction"_a,
             py::keep_alive<1, 3>());
}