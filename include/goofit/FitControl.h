#pragma once

#include "goofit/GlobalCudaDefines.h"

#include <vector>
#include <string>

class PdfBase;

class FitControl {
public:
    FitControl(bool bin, std::string mn);
    ~FitControl();

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
    void setOwner(PdfBase* dat);

protected:
    bool errorsOnBins {false};

private:
    bool binned;
    std::string metricName;
    PdfBase* owner {nullptr};
};

class UnbinnedNllFit : public FitControl {
public:
    UnbinnedNllFit();
};

class BinnedNllFit : public FitControl {
public:
    BinnedNllFit();
};

class BinnedErrorFit : public FitControl {
public:
    BinnedErrorFit();
};

class BinnedChisqFit : public FitControl {
public:
    BinnedChisqFit();
};

