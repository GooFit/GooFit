#pragma once

#include <goofit/Error.h>
#include <goofit/GlobalCudaDefines.h>

#include <string>
#include <vector>

namespace GooFit {

enum class EvalFunc : size_t { Eval = 0, NLL, Prob, BinAvg, BinWithError, Chisq };

static const char *evalfunc_vals[]
    = {"ptr_to_Eval", "ptr_to_NLL", "ptr_to_Prob", "ptr_to_BinAvg", "ptr_to_BinWithError", "ptr_to_Chisq"};

static auto evalfunc_to_string(EvalFunc val) -> const char * { return evalfunc_vals[static_cast<size_t>(val)]; }

class PdfBase;

class FitControl {
  public:
    FitControl(bool bin, EvalFunc mn)
        : binned(bin)
        , metricFunc(mn) {}

    inline auto binnedFit() const -> bool { return binned; }
    inline auto binErrors() const -> bool { return errorsOnBins; }
    inline auto metricIsPdf() const -> bool { return !errorsOnBins; }
    inline auto getMetric() const -> EvalFunc { return metricFunc; }
    inline auto getName() const -> std::string { return evalfunc_to_string(getMetric()); }

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
