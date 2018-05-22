#pragma once

#include <goofit/Error.h>
#include <goofit/GlobalCudaDefines.h>

#include <string>
#include <vector>

namespace GooFit {

enum class EvalFunc : size_t { Eval = 0, NLL, Prob, BinAvg, BinWithError, Chisq };

constexpr const char *evalfunc_vals[]
    = {"ptr_to_Eval", "ptr_to_NLL", "ptr_to_Prob", "ptr_to_BinAvg", "ptr_to_BinWithError", "ptr_to_Chisq"};

constexpr const char *evalfunc_to_string(EvalFunc val) { return evalfunc_vals[static_cast<size_t>(val)]; }

class PdfBase;

class FitControl {
  public:
    FitControl(bool bin, EvalFunc mn)
        : binned(bin)
        , metricFunc(mn) {}

    inline bool binnedFit() const { return binned; }
    inline bool binErrors() const { return errorsOnBins; }
    inline bool metricIsPdf() const { return !errorsOnBins; }
    inline EvalFunc getMetric() const { return metricFunc; }
    inline std::string getName() const { return evalfunc_to_string(getMetric()); }

  protected:
    bool errorsOnBins{false};

  private:
    bool binned;
    EvalFunc metricFunc;
};

class UnbinnedNllFit : public FitControl {
  public:
    UnbinnedNllFit()
        : FitControl(false, EvalFunc::NLL) {}
};

class BinnedNllFit : public FitControl {
  public:
    BinnedNllFit()
        : FitControl(true, EvalFunc::BinAvg) {}
};

class BinnedErrorFit : public FitControl {
  public:
    BinnedErrorFit()
        : FitControl(true, EvalFunc::BinWithError) {
        errorsOnBins = true;
    }
};

class BinnedChisqFit : public FitControl {
  public:
    BinnedChisqFit()
        : FitControl(true, EvalFunc::Chisq) {}
};

class ProbFit : public FitControl {
  public:
    ProbFit()
        : FitControl(false, EvalFunc::Prob) {}
};

class EvalFit : public FitControl {
  public:
    EvalFit()
        : FitControl(false, EvalFunc::Eval) {}
};

} // namespace GooFit
