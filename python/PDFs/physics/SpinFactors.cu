#include <pybind11/pybind11.h>

#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/SpinFactors.h>
#include <goofit/Variable.h>

using namespace GooFit;
namespace py = pybind11;

void init_SpinFactors(py::module &m) {
    py::enum_<SF_4Body>(m, "SF_4Body", py::arithmetic())
            .value("", SF_4Body::DtoPP1_PtoSP2_StoP3P4)
            .value("", SF_4Body::DtoPP1_PtoVP2_VtoP3P4)
            .value("", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_S)
            .value("", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_P)
            .value("", SF_4Body::DtoV1V2_V1toP1P2_V2toP3P4_D)
            .value("", SF_4Body::DtoAP1_AtoVP2_VtoP3P4)
            .value("", SF_4Body::DtoAP1_AtoVP2Dwave_VtoP3P4)
            .value("", SF_4Body::DtoVS_VtoP1P2_StoP3P4)
            .value("", SF_4Body::DtoV1P1_V1toV2P2_V2toP3P4)
            .value("", SF_4Body::DtoAP1_AtoSP2_StoP3P4)
            .value("", SF_4Body::DtoTP1_TtoVP2_VtoP3P4)
            .value("", SF_4Body::FF_12_34_L1)
            .value("", SF_4Body::FF_12_34_L2)
            .value("", SF_4Body::FF_123_4_L1)
            .value("", SF_4Body::FF_123_4_L2)
            .value("", SF_4Body::ONE);

    py::class_<SpinFactor, GooPdf>(m, "SpinFactor")
        .def(py::init<std::string, SF_4Body, unsigned int, unsigned int, unsigned int, unsigned int>())

        ;
}
