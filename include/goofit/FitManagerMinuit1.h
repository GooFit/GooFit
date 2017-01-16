#ifndef FITMANAGER_MINUIT1_HH
#define FITMANAGER_MINUIT1_HH

#include "TMinuit.h"

extern PdfBase* pdfPointer;
extern int numPars;

void FitFun(int& npar, double* gin, double& fun, double* fp, int iflag);

class FitManager {
public:
    FitManager(PdfBase* dat);
    ~FitManager();
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
    void setupMinuit();
    void runMigrad();
    void fit();
    TMinuit* getMinuitObject() {
        return minuit;
    }
    void getMinuitValues() const;
    void getMinuitStatus(double& fmin, double& fedm, double& errdef, int& npari, int& nparx, int& istat) const;
    TMinuit* minuit;
private:
    double overrideCallLimit;
    bool _useHesseBefore;
    bool _useHesse;
    bool _useMinos;
    bool _useImprove;
};

#endif
