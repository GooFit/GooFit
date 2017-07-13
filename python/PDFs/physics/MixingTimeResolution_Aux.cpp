#pragma once
#include <pybind11/pybind11.h>

#include <goofit/PDFs/physics/MixingTimeResolution_Aux.h>

using namespace GooFit;
namespace py = pybind11;

template<class Base = MixingTimeResolution>

class PyMixingTimeResolution : public Base {

public:

    using Base::Base ;

    fptype normalisation(fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing, fptype ymixing) const override {
        PYBIND11_OVERLOAD_PURE(
                    fptype,
                    Base,
                    normalisation,
                    di1, di2, di3, di4, tau, xmixing, ymixing
                    );
    }

    void createParameters(std::vector<unsigned int> &pindices, PdfBase *dis) override {
        PYBIND11_OVERLOAD_PURE(
                    void,
                    Base,
                    createParameters,
                    pindices,dis
                    );
    }

};
void init_MixingTimeResolution(py::module &m) {
    py::class_<MixingTimeResolution, PyMixingTimeResolution<>>(m, "MixingTimeResolution")
            .def(py::init<>())
            .def("createParameters", &MixingTimeResolution::createParameters)
            .def("normalisation", &MixingTimeResolution::normalisation)
    ;
}

