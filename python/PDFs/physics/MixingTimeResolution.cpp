#include <goofit/Python.h>

#include <goofit/PDFs/physics/MixingTimeResolution.h>
#include <goofit/docs/PDFs/physics/MixingTimeResolution.h>

using namespace GooFit;

template <class Base = GooFit::MixingTimeResolution>
class PyMixingTimeResolution : public Base {
  public:
    using Base::Base;

    GooFit::fptype normalization(GooFit::fptype di1,
                                 GooFit::fptype di2,
                                 GooFit::fptype di3,
                                 GooFit::fptype di4,
                                 GooFit::fptype tau,
                                 GooFit::fptype xmixing,
                                 GooFit::fptype ymixing) const override {
        PYBIND11_OVERLOAD_PURE(GooFit::fptype, Base, normalization, di1, di2, di3, di4, tau, xmixing, ymixing);
    }

    void createParameters(GooFit::PdfBase *dis) override { PYBIND11_OVERLOAD_PURE(void, Base, createParameters, dis); }
};

void init_MixingTimeResolution(py::module &m) {
    py::class_<MixingTimeResolution, PyMixingTimeResolution<>>(m, "MixingTimeResolution")
        .def(py::init<std::string>(), "n"_a)

        .def("createParameters", &MixingTimeResolution::createParameters)
        .def("normalization", &MixingTimeResolution::normalization)

        .def_static("help", []() { return HelpPrinter(MixingTimeResolution_docs); })

        ;
}
