#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/physics/LineshapesPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_LineshapesPdf(py::module &m) {

    py::enum_<LS>(m, "LS", py::arithmetic())
            .value("", LS::ONE)
            .value("", LS::BW)
            .value("", LS::Lass)
            .value("", LS::Lass_M3)
            .value("", LS::nonRes)
            .value("", LS::Bugg)
            .value("", LS::Bugg3)
            .value("", LS::Flatte)
            .value("", LS::SBW);

    py::enum_<FF>(m, "FF", py::arithmetic())
            .value("", FF::One)
            .value("", FF::BL)
            .value("", FF::BL_Prime)
            .value("", FF::BL2);

    py::class_<Lineshape, GooPdf>(m, "Lineshape")
        .def(py::init<std::string,
                      Variable *,
                      Variable *,
                      unsigned int,
                      unsigned int,
                      LS,
                      FF,
                      fptype,
                      std::vector<Variable *>>())

        ;
}
