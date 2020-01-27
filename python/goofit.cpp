#include <goofit/Python.h>

#include <pybind11/eval.h>
#include <pybind11/pybind11.h>

void init_HelpPrinter(py::module &);
void init_DataSet(py::module &);
void init_BinnedDataSet(py::module &);
void init_UnbinnedDataSet(py::module &);
void init_Variable(py::module &);
void init_FitManager(py::module &);
void init_PdfBase(py::module &);
void init_GooPdf(py::module &);
void init_CombinePdf(py::module &);
void init_Version(py::module &);
void init_FitControl(py::module &);
void init_Application(py::module &);

// Basic
void init_ArgusPdf(py::module &);
void init_BifurGaussPdf(py::module &);
void init_BernsteinPdf(py::module &);
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
void init_DalitzPlotHelpers(py::module &);
void init_AmpNBodyBase(py::module &);
void init_Amp3BodyBase(py::module &);
void init_Amp3Body(py::module &);
void init_Amp3Body_TD(py::module &);
void init_Amp4BodyBase(py::module &);
void init_Amp4Body(py::module &);
void init_Amp4Body_TD(py::module &);
void init_DalitzPlotter(py::module &);
void init_DalitzVetoPdf(py::module &);
void init_Amp3Body_IS(py::module &);
void init_Lineshapes(py::module &);
void init_MixingTimeResolution(py::module &);
void init_ResonancePdf(py::module &);
void init_SpinFactors(py::module &);
void init_ThreeGaussResolution(py::module &);
void init_TruthResolution(py::module &);
void init_SquareDalitzEffPdf(py::module &);

// Utilities
void init_VariableBinTransform1DPdf(py::module &);

PYBIND11_MODULE(_goofit, m) {
    m.doc() = "Python interface for GooFit";

    py::module::import("goofit.minuit2");

    init_HelpPrinter(m);
    init_Variable(m);
    init_DataSet(m);
    init_BinnedDataSet(m);
    init_UnbinnedDataSet(m);
    init_FitManager(m);
    init_PdfBase(m);
    init_GooPdf(m);
    init_CombinePdf(m);
    init_Version(m);
    init_FitControl(m);
    init_Application(m);

    // Basic
    init_ArgusPdf(m);
    init_BifurGaussPdf(m);
    init_BinTransformPdf(m);
    init_BernsteinPdf(m);
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
    init_DalitzPlotHelpers(m);
    init_DalitzPlotter(m);
    init_AmpNBodyBase(m);
    init_Amp3BodyBase(m);
    init_Amp3Body(m);
    init_Amp3Body_TD(m);
    init_Amp4BodyBase(m);
    init_Amp4Body(m);
    init_Amp4Body_TD(m);
    init_DalitzVetoPdf(m);
    init_Amp3Body_IS(m);
    init_Lineshapes(m);
    init_MixingTimeResolution(m);
    init_ResonancePdf(m);
    init_SpinFactors(m);
    init_ThreeGaussResolution(m);
    init_TruthResolution(m);
    init_SquareDalitzEffPdf(m);

    // Utilities
    init_VariableBinTransform1DPdf(m);

    // Setup for iPython
    py::object ip = py::none();

    try {
        py::object ipython = py::module::import("IPython");
        ip                 = ipython.attr("get_ipython")();
    } catch(const py::error_already_set &) {
    }

    if(!ip.is_none()) {
        py::object html_formatter = ip.attr("display_formatter").attr("formatters")["text/markdown"];
        auto locals               = py::dict("html_formatter"_a = html_formatter, "m"_a = m);
        py::exec("html_formatter.for_type(type(m.PdfBase), lambda x: repr(x.help()) if hasattr(x, 'help') else None)",
                 py::globals(),
                 locals);
    }
}
