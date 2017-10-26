#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/physics/LineshapesPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;
using namespace pybind11::literals;

void init_LineshapesPdf(py::module &m) {
    py::enum_<FF>(m, "FF", py::arithmetic())
        .value("One", FF::One)
        .value("BL", FF::BL)
        .value("BL_Prime", FF::BL_Prime)
        .value("BL2", FF::BL2);

    auto m_ls = m.def_submodule("Lineshapes");

    py::class_<Lineshapes::One>(m_ls, "One")
        .def(py::init<std::string, Variable *, Variable *, unsigned int, unsigned int, FF, fptype>(),
             "Create a constant lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a = FF::BL_Prime,
             "radius"_a  = 1.5);

    py::class_<Lineshapes::LASS>(m_ls, "LASS")
        .def(py::init<std::string, Variable *, Variable *, unsigned int, unsigned int, FF, fptype>(),
             "Create a LASS lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a = FF::BL_Prime,
             "radius"_a  = 1.5);

    py::class_<Lineshapes::NonRes>(m_ls, "NonRes")
        .def(py::init<std::string, Variable *, Variable *, unsigned int, unsigned int, FF, fptype>(),
             "Create a Non Resonant lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a = FF::BL_Prime,
             "radius"_a  = 1.5);

    py::class_<Lineshapes::Flatte>(m_ls, "Flatte")
        .def(py::init<std::string, Variable *, Variable *, unsigned int, unsigned int, FF, fptype>(),
             "Create a Flatte lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a = FF::BL_Prime,
             "radius"_a  = 1.5);

    py::class_<Lineshapes::Bugg>(m_ls, "Bugg")
        .def(py::init<std::string, Variable *, Variable *, unsigned int, unsigned int, FF, fptype>(),
             "Create a constant lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a = FF::BL_Prime,
             "radius"_a  = 1.5);

    py::class_<Lineshapes::Bugg3>(m_ls, "Bugg3")
        .def(py::init<std::string, Variable *, Variable *, unsigned int, unsigned int, FF, fptype>(),
             "Create a constant lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a = FF::BL_Prime,
             "radius"_a  = 1.5);

    py::class_<Lineshapes::RBW>(m_ls, "RBW")
        .def(py::init<std::string, Variable *, Variable *, unsigned int, unsigned int, FF, fptype>(),
             "Create a RBW lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a = FF::BL_Prime,
             "radius"_a  = 1.5);

    py::class_<Lineshapes::SBW>(m_ls, "SBW")
        .def(py::init<std::string, Variable *, Variable *, unsigned int, unsigned int, FF, fptype>(),
             "Create a SBW lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a = FF::BL_Prime,
             "radius"_a  = 1.5);

    py::class_<Lineshapes::GLASS>(m_ls, "GLASS")
        .def(py::init<std::string,
                      Variable *,
                      Variable *,
                      unsigned int,
                      unsigned int,
                      FF,
                      fptype,
                      std::vector<Variable *>>(),
             "Create a G-LASS lineshape",

             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a        = FF::BL_Prime,
             "radius"_a         = 1.5,
             "AdditionalVars"_a = std::vector<Variable *>());

    py::class_<Lineshapes::GSpline>(m_ls, "GSpline")
        .def(py::init<std::string,
                      Variable *,
                      Variable *,
                      unsigned int,
                      unsigned int,
                      FF,
                      fptype,
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
