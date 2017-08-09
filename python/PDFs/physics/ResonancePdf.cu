#include <pybind11/pybind11.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/physics/ResonancePdf.h>

using namespace GooFit;
namespace py = pybind11;
using namespace pybind11::literals;

void init_ResonancePdf(py::module &m) {
    m.attr("MAXNKNOBS") = MAXNKNOBS;
    
    py::enum_<ResPdfType>(m, "ResPdfType", py::arithmetic())
        .value("RBW", ResPdfType::RBW)
        .value("LASS", ResPdfType::LASS)
        .value("GS", ResPdfType::GS)
        .value("FLATTE", ResPdfType::FLATTE)
        .value("GAUSS", ResPdfType::GAUSS)
        .value("SPLINE", ResPdfType::SPLINE)
        .value("NONRES", ResPdfType::NONRES)
        .export_values();
    
    py::class_<ResonancePdf, GooPdf>(m, "ResonancePdf")
        .def(py::init<std::string, ResPdfType, Variable *, Variable *, Variable *, Variable *, unsigned int, unsigned int, bool>(),
             "Constructor for regular BW,Gounaris-Sakurai,LASS",
             "name"_a, "ResPdfType"_a, "ar"_a, "ai"_a, "mass"_a, "width"_a, "sp"_a, "cyc"_a, "symmDP"_a = false)
    
        .def(py::init<std::string, ResPdfType, Variable *, Variable *>(), "Constructor for regular BW,Gounaris-Sakurai,LASS"
            "Constructor for NONRES",
            "name"_a, "ResPdfType"_a, "ar"_a, "ai"_a)
    
        .def(py::init<std::string, ResPdfType, Variable *, Variable *, Variable *, Variable *, unsigned int, bool>(),
             "Constructor for regular GAUSS",
             "name"_a, "ResPdfType"_a, "ar"_a, "ai"_a, "mean"_a, "sigma"_a, "cyc"_a, "symmDP"_a = false)
        
        .def(py::init<std::string, ResPdfType, Variable *, Variable *, Variable *, Variable*, Variable *, unsigned int, bool>(),
             "Constructor for regular FLATTE",
             "name"_a, "ResPdfType"_a, "ar"_a, "ai"_a, "mean"_a, "g1"_a, "rg2og1"_a, "cyc"_a, "symmDP"_a = false)
        
        .def(py::init<std::string, ResPdfType, Variable *, Variable *, std::vector<fptype>&, std::vector<Variable*>&, std::vector<Variable*>&, unsigned int, bool>(),
             "Constructor for regular CUBIC spline",
             "name"_a, "ResPdfType"_a, "ar"_a, "ai"_a, "HH_bin_limits"_a, "pwa_coefs_reals"_a, "pwa_coefs_imags"_a, "cyc"_a, "symmDP"_a = false)
        ;
}








