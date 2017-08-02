#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_DataSet(py::module &);
void init_BinnedDataSet(py::module &);
void init_UnbinnedDataSet(py::module &);
void init_Variable(py::module &);
void init_FitManager(py::module &);
void init_PdfBase(py::module &);
void init_GooPdf(py::module &);
void init_Version(py::module &);

// Basic
void init_ArgusPdf(py::module &);
void init_BifurGaussPdf(py::module &);
void init_BinTransformPdf(py::module &);
void init_BWPdf(py::module &);
void init_CorrGaussianPdf(py::module &);
void init_CrystalBallPdf(py::module &);
void init_ExpGausPdf(py::module &);
void init_ExpPdf(py::module &);
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

// Combine
void init_AddPdf(py::module &);
void init_CompositePdf(py::module &);
void init_ConvolutionPdf(py::module &);
void init_EventWeightedAddPdf(py::module &);
void init_MappedPdf(py::module &);
void init_ProdPdf(py::module &);

// Physics
void init_DalitzPlotPdf(py::module &);
void init_DalitzVetoPdf(py::module &);
void init_DP4Pdf(py::module &);
void init_IncoherentSumPdf(py::module &);
void init_LineshapesPdf(py::module &);
void init_MixingTimeResolution(py::module &);
void init_ResonancePdf(py::module &);
void init_SpinFactors(py::module &);
void init_Tddp4Pdf(py::module &);
void init_TddpPdf(py::module &);
void init_ThreeGaussResolution(py::module &);
void init_TruthResolution(py::module &);

// Utilities


PYBIND11_PLUGIN(_goofit) {
    py::module m("_goofit", "Python interface for GooFit");

    init_Variable(m);
    init_DataSet(m);
    init_BinnedDataSet(m);
    init_UnbinnedDataSet(m);
    init_FitManager(m);
    init_PdfBase(m);
    init_GooPdf(m);
    init_Version(m);

    // Basic
    init_ArgusPdf(m);
    init_BifurGaussPdf(m);
    init_BinTransformPdf(m);
    init_BWPdf(m);
    init_CorrGaussianPdf(m);
    init_CrystalBallPdf(m);
    init_ExpGausPdf(m);
    init_ExpPdf(m);
    init_GaussianPdf(m);
    init_InterHistPdf(m);
    init_JohnsonSUPdf(m);
    init_KinLimitBWPdf(m);
    init_LandauPdf(m);
    init_NovosibirskPdf(m);
    init_PolynomialPdf(m);
    init_ScaledGaussianPdf(m);
    init_SmoothHistogramPdf(m);
    init_StepPdf(m);
    init_TrigThresholdPdf(m);
    init_VoigtianPdf(m);
    
    // Combine
    init_AddPdf(m);
    init_CompositePdf(m);
    init_ConvolutionPdf(m);
    init_EventWeightedAddPdf(m);
    init_MappedPdf(m);
    init_ProdPdf(m);
    
    // Physics
    init_DalitzPlotPdf(m);
    init_DalitzVetoPdf(m);
    init_DP4Pdf(m);
    init_IncoherentSumPdf(m);
    init_LineshapesPdf(m);
    init_MixingTimeResolution(m);
    init_ResonancePdf(m);
    init_SpinFactors(m);
    init_Tddp4Pdf(m);
    init_TddpPdf(m);
    init_ThreeGaussResolution(m);
    init_TruthResolution(m);

    // Utilities


    return m.ptr();
}
