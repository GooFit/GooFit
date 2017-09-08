#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/PDFs/physics/LineshapesPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

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
