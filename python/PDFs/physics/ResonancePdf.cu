#include <pybind11/pybind11.h>

#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;
using namespace pybind11::literals;

void init_ResonancePdf(py::module &m) {
    m.attr("MAXNKNOBS") = MAXNKNOBS;

    auto m_ls = m.def_submodule("Lineshapes");
    
    py::class_<ResonancePdf, GooPdf>(m, "ResonancePdf");
    
    py::class_<Resonances::RBW, ResonancePdf>(m_ls, "RBW")
        .def(py::init<std::string,
                      Variable *,
                      Variable *,
                      Variable *,
                      Variable *,
                      unsigned int,
                      unsigned int,
                      bool>(),
             "Constructor for regular BW",
             "name"_a,
             "ar"_a,
             "ai"_a,
             "mass"_a,
             "width"_a,
             "sp"_a,
             "cyc"_a,
             "symmDP"_a = false)
    ;
    
    py::class_<Resonances::GS, ResonancePdf>(m, "GS")
    .def(py::init<std::string,
         Variable *,
         Variable *,
         Variable *,
         Variable *,
         unsigned int,
         unsigned int,
         bool>(),
         "Constructor for regular Gounaris-Sakurai",
         "name"_a,
         "ar"_a,
         "ai"_a,
         "mass"_a,
         "width"_a,
         "sp"_a,
         "cyc"_a,
         "symmDP"_a = false)
    ;

    py::class_<Resonances::LASS, ResonancePdf>(m, "LASS")
    .def(py::init<std::string,
         Variable *,
         Variable *,
         Variable *,
         Variable *,
         unsigned int,
         unsigned int,
         bool>(),
         "Constructor for LASS",
         "name"_a,
         "ar"_a,
         "ai"_a,
         "mass"_a,
         "width"_a,
         "sp"_a,
         "cyc"_a,
         "symmDP"_a = false)
    ;
    
    py::class_<Resonances::NonRes, ResonancePdf>(m, "NonRes")
        .def(py::init<std::string, Variable *, Variable *>(),
             "Constructor for NonResonant",
             "name"_a,
             "ar"_a,
             "ai"_a);

    py::class_<Resonances::Gauss, ResonancePdf>(m, "Gauss")
        .def(py::init<std::string, Variable *, Variable *, Variable *, Variable *, unsigned int>(),
             "Constructor for regular GAUSS",
             "name"_a,
             "ar"_a,
             "ai"_a,
             "mean"_a,
             "sigma"_a,
             "cyc"_a);

    py::class_<Resonances::FLATTE, ResonancePdf>(m, "FLATTE")
        .def(py::init<std::string,
                      Variable *,
                      Variable *,
                      Variable *,
                      Variable *,
                      Variable *,
                      unsigned int,
                      bool>(),
             "Constructor for regular FLATTE",
             "name"_a,
             "ar"_a,
             "ai"_a,
             "mean"_a,
             "g1"_a,
             "rg2og1"_a,
             "cyc"_a,
             "symmDP"_a);

    py::class_<Resonances::Spline, ResonancePdf>(m, "Spline")
        .def(py::init<std::string,
                      Variable *,
                      Variable *,
                      std::vector<fptype> &,
                      std::vector<Variable *> &,
                      std::vector<Variable *> &,
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
             "symmDP"_a = false);
}
