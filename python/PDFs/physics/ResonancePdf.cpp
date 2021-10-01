#include <goofit/Python.h>

#include <pybind11/stl.h>

#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/Variable.h>
#include <goofit/docs/PDFs/physics/resonances/Resonance.h>

#if GOOFIT_KMATRIX
#include <goofit/PDFs/physics/resonances/kMatrix.h>
#endif

using namespace GooFit;

void init_ResonancePdf(py::module &m) {
    auto m_ls = m.def_submodule("Resonances");

    py::class_<ResonancePdf, GooPdf>(m, "ResonancePdf")
        .def("get_amp_real", &ResonancePdf::get_amp_real)
        .def("get_amp_img", &ResonancePdf::get_amp_img)
        .def_static("help", []() { return HelpPrinter(ResonancePdf_docs); });

    py::class_<Resonances::RBW, ResonancePdf>(m_ls, "RBW")
        .def(py::init<std::string, Variable, Variable, Variable, Variable, unsigned int, unsigned int, bool, bool>(),
             "Constructor for regular BW",
             "name"_a,
             "ar"_a,
             "ai"_a,
             "mass"_a,
             "width"_a,
             "sp"_a,
             "cyc"_a,
             "norm"_a = true,
             "sym"_a  = false);

    py::class_<Resonances::GS, ResonancePdf>(m_ls, "GS")
        .def(py::init<std::string, Variable, Variable, Variable, Variable, unsigned int, unsigned int, bool, bool>(),
             "Constructor for regular Gounaris-Sakurai",
             "name"_a,
             "ar"_a,
             "ai"_a,
             "mass"_a,
             "width"_a,
             "sp"_a,
             "cyc"_a,
             "norm"_a = true,
             "sym"_a  = false);

    py::class_<Resonances::RhoOmegaMix, ResonancePdf>(m_ls, "RhoOmegaMix")
        .def(py::init<std::string,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      unsigned int,
                      unsigned int,
                      bool,
                      bool>(),
             "Constructor for rho omega mixing",
             "name"_a,
             "ar"_a,
             "ai"_a,
             "omega_mass"_a,
             "omega_width"_a,
             "rho_mass"_a,
             "rho_width"_a,
             "real"_a,
             "imag"_a,
             "delta"_a,
             "sp"_a,
             "cyc"_a,
             "norm"_a = true,
             "sym"_a  = false);

    py::class_<Resonances::LASS, ResonancePdf>(m_ls, "LASS")
        .def(py::init<std::string,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      unsigned int,
                      unsigned int,
                      bool>(),
             "Constructor for LASS",
             "name"_a,
             "ar"_a,
             "ai"_a,
             "mass"_a,
             "width"_a,
             "_a"_a,
             "_r"_a,
             "_R"_a,
             "_phiR"_a,
             "_B"_a,
             "_phiB"_a,
             "sp"_a,
             "cyc"_a,
             "norm"_a = true);

#if GOOFIT_KMATRIX
    py::class_<Resonances::kMatrix, ResonancePdf>(m_ls, "kMatrix")
        .def(py::init<std::string,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      Variable,
                      std::vector<Variable>,
                      std::vector<Variable>,
                      std::vector<Variable>,
                      std::vector<Variable>,
                      std::vector<Variable>,
                      std::vector<Variable>,
                      unsigned int,
                      unsigned int>(),
             "kMatrix resonance pdf",
             "name"_a,
             "a_r"_a,
             "a_i"_a,
             "sA0"_a,
             "sA"_a,
             "s0_prod"_a,
             "s0_scatt"_a,
             "beta_r"_a,
             "beta_i"_a,
             "f_prod_r"_a,
             "f_prod_i"_a,
             "fscat"_a,
             "poles"_a,
             "L"_a,
             "Mpair"_a);
#endif

    py::class_<Resonances::NonRes, ResonancePdf>(m_ls, "NonRes")
        .def(py::init<std::string, Variable, Variable>(), "Constructor for NonResonant", "name"_a, "ar"_a, "ai"_a);

    py::class_<Resonances::Gauss, ResonancePdf>(m_ls, "Gauss")
        .def(py::init<std::string, Variable, Variable, Variable, Variable, unsigned int>(),
             "Constructor for regular GAUSS",
             "name"_a,
             "ar"_a,
             "ai"_a,
             "mean"_a,
             "sigma"_a,
             "cyc"_a);

    py::class_<Resonances::FLATTE, ResonancePdf>(m_ls, "FLATTE")
        .def(py::init<std::string, Variable, Variable, Variable, Variable, Variable, unsigned int, bool>(),
             "Constructor for regular FLATTE",
             "name"_a,
             "ar"_a,
             "ai"_a,
             "mean"_a,
             "g1"_a,
             "rg2og1"_a,
             "cyc"_a,
             "symmDP"_a);

    py::class_<Resonances::Spline, ResonancePdf>(m_ls, "Spline")
        .def(py::init<std::string,
                      Variable,
                      Variable,
                      std::vector<fptype> &,
                      std::vector<Variable> &,
                      std::vector<Variable> &,
                      unsigned int,
                      bool>(),
             "Constructor for regular cubic spline",
             "name"_a,
             "ar"_a,
             "ai"_a,
             "HH_bin_limits"_a,
             "pwa_coefs_reals"_a,
             "pwa_coefs_imags"_a,
             "cyc"_a,
             "symmDP"_a = false,
             py::keep_alive<1, 5>(),
             py::keep_alive<1, 6>(),
             py::keep_alive<1, 7>());
}
