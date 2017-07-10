#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <goofit/Variable.h>
#include <goofit/PDFs/physics/ThreeGaussResolution_Aux.h>
#include "MixingTimeResolution_Aux.h"


using namespace GooFit;
namespace py = pybind11;

template<class ThreeGaussResolutionBase = ThreeGaussResolution> class PyThreeGaussResolution :
        public PyMixingTimeResolution<ThreeGaussResolutionBase> {

    using PyMixingTimeResolution<ThreeGaussResolutionBase>::PyMixingTimeResolution;

    fptype normalisation(fptype di1, fptype di2, fptype di3, fptype di4, fptype tau, fptype xmixing, fptype ymixing) const override {
        PYBIND11_OVERLOAD(
                    fptype,
                    ThreeGaussResolutionBase,
                    normalisation,
                    di1, di2, di3, di4, tau, xmixing, ymixing
                    );
    }

    void createParameters(std::vector<unsigned int> &pindices, PdfBase *dis) override {
        PYBIND11_OVERLOAD(
                    void,
                    ThreeGaussResolutionBase,
                    createParameters,
                    pindices,dis
                    );
    }

};


void init_ThreeGaussResolution(py::module &m) {
    py::class_<ThreeGaussResolution, PyThreeGaussResolution<>>(m, "ThreeGaussResolution")
        .def(py::init<Variable *, Variable *,Variable *, Variable *,Variable *, Variable *,Variable *, Variable *>())
        ;
}








