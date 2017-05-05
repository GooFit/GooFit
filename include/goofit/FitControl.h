#pragma once

#include "goofit/GlobalCudaDefines.h"
#include "goofit/Error.h"

#include <vector>
#include <string>

class PdfBase;

class FitControl {
public:
    FitControl(bool bin, std::string mn)
      : binned(bin)
      , metricName(mn) {}

    inline bool binnedFit() const {
        return binned;
    }
    inline bool binErrors() const {
        return errorsOnBins;
    }
    inline bool metricIsPdf() const {
        return !errorsOnBins;
    }
    inline std::string getMetric() const {
        return metricName;
    }
    inline PdfBase* getOwner() const {
        return owner;
    }
    void setOwner(PdfBase* dat) {
        if(owner == nullptr)
            throw GooFit::GeneralError("Owner is null");
        owner = dat;
    }

protected:
    bool errorsOnBins {false};

private:
    bool binned;
    std::string metricName;
    PdfBase* owner {nullptr};
};

class UnbinnedNllFit : public FitControl {
public:
    UnbinnedNllFit() : FitControl(false, "ptr_to_NLL") {}
};

class BinnedNllFit : public FitControl {
public:
    BinnedNllFit() : FitControl(true, "ptr_to_BinAvg") {}
};

class BinnedErrorFit : public FitControl {
public:
    BinnedErrorFit() : FitControl(true, "ptr_to_BinWithError") {
        errorsOnBins = true;
    }
};

class BinnedChisqFit : public FitControl {
public:
    BinnedChisqFit() : FitControl(true, "ptr_to_Chisq") {}
};

