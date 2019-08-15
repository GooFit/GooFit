#include <goofit/Python.h>

#include <pybind11/stl.h>

#include <goofit/PDFs/physics/Lineshapes.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/Variable.h>

using namespace GooFit;

void init_Lineshapes(py::module &m) {
    py::enum_<FF>(m, "FF", py::arithmetic())
        .value("One", FF::One)
        .value("BL", FF::BL)
        .value("BL_Prime", FF::BL_Prime)
        .value("BL2", FF::BL2)

        ;

    auto m_ls = m.def_submodule("Lineshapes");

    py::class_<Lineshape, GooPdf>(m, "Lineshape");

    py::class_<Lineshapes::One, Lineshape>(m_ls, "One")
        .def(py::init<std::string, Variable, Variable, unsigned int, unsigned int, FF, fptype>(),
             "Create a constant lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a = FF::BL_Prime,
             "radius"_a  = 1.5);

    py::class_<Lineshapes::LASS, Lineshape>(m_ls, "LASS")
        .def(py::init<std::string, Variable, Variable, unsigned int, unsigned int, FF, fptype>(),
             "Create a LASS lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a = FF::BL_Prime,
             "radius"_a  = 1.5);

    py::class_<Lineshapes::kMatrix, Lineshape>(m_ls, "kMatrix")
<<<<<<< HEAD
      .def(py::init<std::string,
                    unsigned int,
                    bool,
                    Variable,
                    Variable,
                    Variable,
                    Variable,
	            std::array<Variable, 5>,
	            std::array<Variable, 5*6>,
	            Variable, 
	            Variable,
                    unsigned int,
                    unsigned int,
                    FF,
                    fptype>(),
           "Create a kMatrix lineshape",
           "name"_a,
           "pterm"_a,
           "is_pole"_a,
           "sA0"_a,
           "sA"_a,
           "s0_prod"_a,
           "s0_scatt"_a,
           "fscat"_a,
           "poles"_a,
           "mass"_a,
           "width"_a,
           "L"_a,
           "Mpair"_a,
           "FormFac"_a = FF::BL_Prime,
           "radius"_a  = 1.5);
=======
        .def(py::init<std::string,
                      unsigned int,
                      bool,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      std::array<Variable, 5>,
                      std::array<Variable, 5 * 6>,
                      Variable,
                      Variable,
                      unsigned int,
                      unsigned int,
                      FF,
                      fptype>(),
             "Create a kMatrix lineshape",
             "name"_a,
             "pterm"_a,
             "is_pole"_a,
             "sA0"_a,
             "sA"_a,
             "s0_prod"_a,
             "s0_scatt"_a,
             "fscat"_a,
             "poles"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a = FF::BL_Prime,
             "radius"_a  = 1.5);
>>>>>>> 19023a91376e86204a3421484ed440ba6019eda5

    py::class_<Lineshapes::NonRes, Lineshape>(m_ls, "NonRes")
        .def(py::init<std::string, Variable, Variable, unsigned int, unsigned int, FF, fptype>(),
             "Create a Non Resonant lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a = FF::BL_Prime,
             "radius"_a  = 1.5);

    py::class_<Lineshapes::Flatte, Lineshape>(m_ls, "Flatte")
        .def(py::init<std::string, Variable, Variable, unsigned int, unsigned int, FF, fptype>(),
             "Create a Flatte lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a = FF::BL_Prime,
             "radius"_a  = 1.5);

    py::class_<Lineshapes::Bugg, Lineshape>(m_ls, "Bugg")
        .def(py::init<std::string, Variable, Variable, unsigned int, unsigned int, FF, fptype>(),
             "Create a Bugg lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a = FF::BL_Prime,
             "radius"_a  = 1.5);

    py::class_<Lineshapes::Bugg3, Lineshape>(m_ls, "Bugg3")
        .def(py::init<std::string, Variable, Variable, unsigned int, unsigned int, FF, fptype>(),
             "Create a Bugg3 lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a = FF::BL_Prime,
             "radius"_a  = 1.5);

    py::class_<Lineshapes::RBW, Lineshape>(m_ls, "RBW")
        .def(py::init<std::string, Variable, Variable, unsigned int, unsigned int, FF, fptype>(),
             "Create a RBW lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a = FF::BL_Prime,
             "radius"_a  = 1.5);

    py::class_<Lineshapes::SBW, Lineshape>(m_ls, "SBW")
        .def(py::init<std::string, Variable, Variable, unsigned int, unsigned int, FF, fptype>(),
             "Create a SBW lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a = FF::BL_Prime,
             "radius"_a  = 1.5);

    py::class_<Lineshapes::GLASS, Lineshape>(m_ls, "GLASS")
        .def(py::init<std::string, Variable, Variable, unsigned int, unsigned int, FF, fptype, std::vector<Variable>>(),
             "Create a G-LASS lineshape",

             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a        = FF::BL_Prime,
             "radius"_a         = 1.5,
             "AdditionalVars"_a = std::vector<Variable>());

    py::class_<Lineshapes::GSpline, Lineshape>(m_ls, "GSpline")
        .def(py::init<std::string,
                      Variable,
                      Variable,
                      unsigned int,
                      unsigned int,
                      FF,
                      fptype,
                      std::vector<Variable>,
                      Lineshapes::spline_t>(),
             "Create a GSpline lineshape",
             "name"_a,
             "mass"_a,
             "width"_a,
             "L"_a,
             "Mpair"_a,
             "FormFac"_a        = FF::BL_Prime,
             "radius"_a         = 1.5,
             "AdditionalVars"_a = std::vector<Variable>(),
             "SplineInfo"_a     = Lineshapes::spline_t(0.0, 0.0, 0));

    py::class_<Amplitude>(m, "Amplitude")
        .def(py::init<std::string,
                      Variable,
                      Variable,
                      std::vector<Lineshape *>,
                      std::vector<SpinFactor *>,
                      unsigned int>(),
             "uniqueDecayStr"_a,
             "ar"_a,
             "ai"_a,
             "LS"_a,
             "SF"_a,
             "nPerm"_a = 1,
             py::keep_alive<1, 5>(),
             py::keep_alive<1, 6>());
}
