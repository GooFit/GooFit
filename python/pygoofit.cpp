#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_BinnedDataSet(py::module &);
void init_Variable(py::module &);

PYBIND11_PLUGIN(pygoofit) {
        py::module m("pygoofit", "Python interface for GooFit");

        init_Variable(m);
        init_BinnedDataSet(m);

        return m.ptr();
}


