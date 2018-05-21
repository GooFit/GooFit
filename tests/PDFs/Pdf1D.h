#pragma once

#include <goofit/FitManager.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Variable.h>

#include <goofit/Catch.h>

#include <memory>

class PdfTest1D {
  protected:
    GooFit::Observable xvar{"xvar", 0, 1};
    std::unique_ptr<GooFit::GooPdf> pdf;
    GooFit::UnbinnedDataSet data{xvar};
    bool result_;

  public:
    PdfTest1D() {}

    void fill() {
        pdf->setData(&data);
        pdf->fillMCDataSimple(100000, /* seed */ 42);
    }

    std::string fit() {
        // Fit to data
        GooFit::FitManagerMinuit2 fitter{pdf.get()};
        fitter.setVerbosity(2);
        std::string fitter_output = capture_stdout([&fitter]() { fitter.fit(); });
        result_                   = fitter;
        return fitter_output;
    }

    bool result() const { return result_; }

    template <typename T, typename... Args>
    void set_pdf(Args... args) {
        pdf.reset(new T("pdf", xvar, args...));
    }
};
