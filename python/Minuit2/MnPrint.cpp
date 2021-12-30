#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <Minuit2/MnPrint.h>

#include <goofit/Version.h>

#ifdef MATHCORE_STANDALONE
#define ROOT_VERSION(x, y, z) 0
#else
#include <RVersion.h>
#endif

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ROOT::Minuit2;

void init_MnPrint(py::module &m) {
#if !defined(MATHCORE_STANDALONE) && GOOFIT_ROOT_FOUND && ROOT_VERSION_CODE < ROOT_VERSION(6, 24, 0)
    auto m_ls = m.def_submodule("MnPrint");

    m_ls.def("GlobalLevel", &MnPrint::Level);
    m_ls.def("SetGlobalLevel", &MnPrint::SetLevel);

#else
    py::class_<MnPrint> m_ls(m, "m_ls");

    m_ls.def_static("GlobalLevel", &MnPrint::GlobalLevel);
    m_ls.def_static("SetGlobalLevel", &MnPrint::SetGlobalLevel);
#endif
}
