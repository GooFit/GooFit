#include <pybind11/pybind11.h>
#include <goofit/PDFs/physics/MixingTimeResolution_Aux.h>


using namespace GooFit;
namespace py = pybind11;

template<class Base = GooFit::MixingTimeResolution>
class PyMixingTimeResolution : public Base {

public:

    using Base::Base ;

    GooFit::fptype normalisation(GooFit::fptype di1, GooFit::fptype di2, GooFit::fptype di3, GooFit::fptype di4, GooFit::fptype tau, GooFit::fptype xmixing, GooFit::fptype ymixing) const override {
        PYBIND11_OVERLOAD_PURE(
                    GooFit::fptype,
                    Base,
                    normalisation,
                    di1, di2, di3, di4, tau, xmixing, ymixing
                    );
    }

    void createParameters(std::vector<unsigned int> &pindices, GooFit::PdfBase *dis) override {
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

            .def("createParameters", &MixingTimeResolution::createParameters)
            .def("normalisation", &MixingTimeResolution::normalisation)
            ;
}

