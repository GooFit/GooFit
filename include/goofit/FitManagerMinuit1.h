#pragma once

#include <memory>
#include <TMinuit.h>

class PdfBase;
extern PdfBase* pdfPointer;

void FitFun(int& npar, double* gin, double& fun, double* fp, int iflag);

namespace GooFit {

class FitManagerMinuit1 {
public:
    FitManagerMinuit1(PdfBase* dat) {
        pdfPointer = dat;
    }
        
    void setMaxCalls(double mxc) {
        overrideCallLimit = mxc;
    }
    void useHesseBefore(bool use=true) {
        _useHesseBefore = use;
    }
    void useHesse(bool use=true) {
        _useHesse = use;
    }
    void useMinos(bool use=true) {
        _useMinos = use;
    }
    void useImprove(bool use=true) {
        _useImprove = use;
    }
    
    // This is expected to prepare the fitting
    void setup();
        
    void runMigrad();
        
    void runFit();
    
    // This runs the fit
    void fit() {
        setup();
        runFit();
    };
        
    TMinuit* getMinuitObject() {
        return minuit_.get();
    }
    void getMinuitValues() const;
        
    void getMinuitStatus(double& fmin, double& fedm, double& errdef, int& npari, int& nparx, int& istat) const;
    
        
private:
    double overrideCallLimit {-1};
    bool _useHesseBefore {true};
    bool _useHesse {true};
    bool _useMinos {false};
    bool _useImprove {false};
    std::unique_ptr<TMinuit> minuit_;
};

}
