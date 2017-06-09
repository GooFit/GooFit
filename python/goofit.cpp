#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_DataSet(py::module &);
void init_BinnedDataSet(py::module &);
void init_UnbinnedDataSet(py::module &);
void init_Variable(py::module &);
void init_PdfBase(py::module &);
void init_GooPdf(py::module &);
void init_ExpPdf(py::module &);
void init_FitManager(py::module &);

PYBIND11_PLUGIN(goofit) {
    py::module m("goofit", "Python interface for GooFit");

    init_Variable(m);
    init_DataSet(m);
    init_BinnedDataSet(m);
    init_UnbinnedDataSet(m);
    init_PdfBase(m);
    init_GooPdf(m);
    init_ExpPdf(m);
    init_FitManager(m);

    return m.ptr();
}
