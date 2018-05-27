#include <goofit/Python.h>

#include <pybind11/stl.h>

#include <goofit/PDFs/basic/BinTransformPdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_BinTransformPdf(py::module &m) {
    py::class_<BinTransformPdf, GooPdf>(m, "BinTransformPdf")
        .def(py::init<std::string,
                      std::vector<Observable>,
                      std::vector<fptype>,
                      std::vector<fptype>,
                      std::vector<int>>(),
             "Returns the global bin of its argument in one dimension.",
             "name"_a,
             "obses"_a,
             "limits"_a,
             "binSizes"_a,
             "numBins"_a);
}
