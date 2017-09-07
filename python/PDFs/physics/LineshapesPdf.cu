#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/physics/LineshapesPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

using namespace pybind11::literals;

void init_LineshapesPdf(py::module &m) {
    py::class_<Lineshape, GooPdf>(m, "Lineshape")
        .def(py::init<std::string,
                      Variable *,
                      Variable *,
                      unsigned int,
                      unsigned int,
                      LS,
                      FF,
                      fptype,
                      std::vector<Variable *>,
                      std::vector<Variable *>,
                      Lineshape::spline_t>(),
             "Create a lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "kind"_a,
             "FormFac"_a,
             "radius"_a,
             "AdditionalVars"_a = std::vector<Variable *>(),
             "Curvatures"_a     = std::vector<Variable *>(),
             "SplineInfo"_a     = Lineshape::spline_t(0.0, 0.0, 0));
}
