#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Minuit2/MnPrint.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ROOT::Minuit2;

void init_MnPrint(py::module &m) {
    auto m_ls = m.def_submodule("MnPrint");

    m_ls.def("Level", &MnPrint::Level);
    m_ls.def("SetLevel", &MnPrint::SetLevel);
}
