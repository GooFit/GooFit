#include <goofit/Python.h>

#include <goofit/PDFs/CombinePdf.h>

using namespace GooFit;

void init_CombinePdf(py::module &m) {
    py::class_<CombinePdf, GooPdf>(m, "CombinePdf")

        ;
}
