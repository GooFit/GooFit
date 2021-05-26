#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Minuit2/MnPrint.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ROOT::Minuit2;

void init_MnPrint(py::module &m) {
    py::class_<MnPrint> m_ls(m, "m_ls");

    m_ls.def_static("GlobalLevel", &MnPrint::GlobalLevel);
    m_ls.def_static("SetGlobalLevel", &MnPrint::SetGlobalLevel);
}
