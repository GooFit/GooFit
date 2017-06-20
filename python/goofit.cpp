#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_DataSet(py::module &);
void init_BinnedDataSet(py::module &);
void init_UnbinnedDataSet(py::module &);
void init_Variable(py::module &);
void init_PdfBase(py::module &);
void init_GooPdf(py::module &);
void init_ArgusPdf(py::module &);
void init_BifurGaussPdf(py::module &);
void init_BinTransformPdf(py::module &);
void init_BWPdf(py::module &);
void init_CorrGaussianPdf(py::module &);
void init_CrystalBallPdf(py::module &);
void init_ExpGausPdf(py::module &);
void init_GaussianPdf(py::module &);
void init_InterHistPdf(py::module &);
void init_JohnsonSUPdf(py::module &);
void init_KinLimitBWPdf(py::module &);
void init_LandauPdf(py::module &);
void init_NovosibirskPdf(py::module &);
void init_PolynomialPdf(py::module &);
void init_ScaledGaussianPdf(py::module &);
void init_SmoothHistogramPdf(py::module &);
void init_StepPdf(py::module &);
void init_TrigThresholdPdf(py::module &);
void init_VoigtianPdf(py::module &);
void init_FitManager(py::module &);

PYBIND11_PLUGIN(goofit) {
    py::module m("goofit", "Python interface for GooFit");

    init_Variable(m);
    init_DataSet(m);
    init_BinnedDataSet(m);
    init_UnbinnedDataSet(m);
    init_PdfBase(m);
    init_GooPdf(m);
    init_ExpPdf(m);
    init_FitManager(m);

    return m.ptr();
}
