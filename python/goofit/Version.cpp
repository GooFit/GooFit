#include <goofit/Python.h>

#include <goofit/Version.h>
#include <goofit/VersionGit.h>

void init_Version(py::module &m) {
    m.attr("__version__")          = py::make_tuple(GOOFIT_VERSION_MAJOR, GOOFIT_VERSION_MINOR, GOOFIT_VERSION_PATCH);
    m.attr("GOOFIT_VERSION_TUPLE") = py::make_tuple(GOOFIT_VERSION_MAJOR, GOOFIT_VERSION_MINOR, GOOFIT_VERSION_PATCH);
    m.attr("GOOFIT_VERSION")       = GOOFIT_VERSION;
    m.attr("GOOFIT_GIT_VERSION")   = GOOFIT_GIT_VERSION;
    m.attr("GOOFIT_TAG")           = GOOFIT_TAG;
    m.attr("GOOFIT_SOURCE_DIR")    = GOOFIT_SOURCE_DIR;

    m.attr("CMAKE_CXX_FLAGS")  = CMAKE_CXX_FLAGS;
    m.attr("CMAKE_BUILD_TYPE") = CMAKE_BUILD_TYPE;
    m.attr("GOOFIT_DEVICE")    = GOOFIT_DEVICE;
    m.attr("GOOFIT_HOST")      = GOOFIT_HOST;
    m.attr("ARCH_FLAGS")       = GOOFIT_ARCH_FLAGS;
    m.attr("ROOT_FOUND")       = GOOFIT_ROOT_FOUND;

    m.attr("INT_MAX") = INT_MAX;
}
