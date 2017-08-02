#include <pybind11/pybind11.h>
#include <goofit/Version.h>

namespace py = pybind11;

void init_Version(py::module &m) {
    m.attr("__version__") = GOOFIT_VERSION;
    m.attr("GOOFIT_MAXPAR") = GOOFIT_MAXPAR;

    m.attr("CMAKE_CXX_FLAGS")   = CMAKE_CXX_FLAGS;
    m.attr("CMAKE_BUILD_TYPE")  = CMAKE_BUILD_TYPE;
    m.attr("GOOFIT_DEVICE")     = GOOFIT_DEVICE;
    m.attr("GOOFIT_HOST")       = GOOFIT_HOST;
    m.attr("ARCH_FLAGS")        = ARCH_FLAGS;
    m.attr("ROOT_FOUND")        = ROOT_FOUND;
}
