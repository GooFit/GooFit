#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/physics/LineshapesPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;
using namespace pybind11::literals;

using namespace pybind11::literals;

void init_LineshapesPdf(py::module &m) {
    py::enum_<LS>(m, "LS", py::arithmetic())
        .value("ONE", LS::ONE)
        .value("BW", LS::BW)
        .value("Lass", LS::Lass)
        .value("Lass_M3", LS::Lass_M3)
        .value("nonRes", LS::nonRes)
        .value("Bugg", LS::Bugg)
        .value("Bugg3", LS::Bugg3)
        .value("Flatte", LS::Flatte)
        .value("SBW", LS::SBW);

    py::enum_<FF>(m, "FF", py::arithmetic())
        .value("One", FF::One)
        .value("BL", FF::BL)
        .value("BL_Prime", FF::BL_Prime)
        .value("BL2", FF::BL2);

    py::class_<Lineshape, GooPdf>(m, "Lineshape")
        .def(py::init<std::string, Variable *, Variable *, unsigned int, unsigned int, LS, FF, fptype>(),
             "Create a lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "kind"_a,
             "FormFac"_a = FF::BL_Prime,
             "radius"_a  = 1.5);

    auto m_ls = m.def_submodule("Lineshapes");

    py::class_<Lineshapes::RBW, Lineshape>(m_ls, "RBW")
        .def(py::init<std::string, Variable *, Variable *, unsigned int, unsigned int, FF, fptype>(),
             "Create a RBW lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a = FF::BL_Prime,
             "radius"_a  = 1.5);

    py::class_<Lineshapes::LASS, Lineshape>(m_ls, "LASS")
        .def(py::init<std::string,
                      Variable *,
                      Variable *,
                      unsigned int,
                      unsigned int,
                      FF,
                      fptype,
                      std::vector<Variable *>>(),
             "Create a LASS lineshape",

             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a        = FF::BL_Prime,
             "radius"_a         = 1.5,
             "AdditionalVars"_a = std::vector<Variable *>());

    py::class_<Lineshapes::GSpline, Lineshape>(m_ls, "GSpline")
        .def(py::init<std::string,
                      Variable *,
                      Variable *,
                      unsigned int,
                      unsigned int,
                      FF,
                      fptype,
                      std::vector<Variable *>,
                      std::vector<Variable *>,
                      Lineshapes::spline_t>(),
             "Create a LASS lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a        = FF::BL_Prime,
             "radius"_a         = 1.5,
             "AdditionalVars"_a = std::vector<Variable *>(),
             "Curvatures"_a     = std::vector<Variable *>(),
             "SplineInfo"_a     = Lineshapes::spline_t(0.0, 0.0, 0));

    py::class_<Amplitude>(m, "Amplitude")
        .def(py::init<std::string,
                      Variable *,
                      Variable *,
                      std::vector<Lineshape *>,
                      std::vector<SpinFactor *>,
                      unsigned int>(),
             "uniqueDecayStr"_a,
             "ar"_a,
             "ai"_a,
             "LS"_a,
             "SF"_a,
             "nPerm"_a = 1,
             py::keep_alive<1, 3>(),
             py::keep_alive<1, 4>());
}
