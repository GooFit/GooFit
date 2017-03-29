#ifndef INCOHERENT_SUM_PDF_HH
#define INCOHERENT_SUM_PDF_HH

#include "goofit/PDFs/GooPdf.h"
#include "goofit/PDFs/TddpPdf.h"
#include "goofit/PDFs/devcomplex.h"

// Very similar class to TddpPdf, but without time dependence
// (so no time resolution or mixing) and ignoring interference between
// waves. This makes the code just different enough, the assumptions are
// just enough changed, that it's not worth trying to modify or subclass
// TddpPdf to deal with both cases. So instead we have a separate
// class with fairly similar structure.


class SpecialIncoherentIntegrator;
class SpecialIncoherentResonanceCalculator;

class IncoherentSumPdf : public GooPdf {
public:
    IncoherentSumPdf(std::string n, Variable* m12, Variable* m13, CountingVariable* eventNumber, DecayInfo* decay,
                     GooPdf* eff);
    // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the
    // incoherent sum. The caching method requires that it be done this way or the ProdPdf
    // normalisation will get *really* confused and give wrong answers.
    __host__ virtual fptype normalise() const;
    __host__ void setDataSize(unsigned int dataSize, unsigned int evtSize = 3);
    __host__ void setForceIntegrals(bool f = true) {
        forceRedoIntegrals = f;
    }

protected:

private:
    DecayInfo* decayInfo;
    Variable* _m12;
    Variable* _m13;
    fptype* dalitzNormRange;

    // Following variables are useful if masses and widths, involved in difficult BW calculation,
    // change infrequently while amplitudes, only used in adding BW results together, change rapidly.
    DEVICE_VECTOR<devcomplex<fptype>>* cachedResonances; // BW (and other) results for each event.
    double* integrals; // Integrals of each BW resonance across the Daliz plot.

    bool* redoIntegral;
    mutable bool forceRedoIntegrals;
    fptype* cachedMasses;
    fptype* cachedWidths;
    int totalEventSize;
    int cacheToUse;
    PdfBase* efficiency;
    SpecialIncoherentIntegrator** integrators;
    SpecialIncoherentResonanceCalculator** calculators;
};

class SpecialIncoherentIntegrator : public thrust::unary_function<thrust::tuple<int, fptype*>, fptype > {
public:
    SpecialIncoherentIntegrator(int pIdx, unsigned int ri);
    EXEC_TARGET fptype operator()(thrust::tuple<int, fptype*> t) const;

private:
    unsigned int resonance_i;
    unsigned int parameters;
};

class SpecialIncoherentResonanceCalculator : public
    thrust::unary_function<thrust::tuple<int, fptype*, int>, devcomplex<fptype>> {
public:

    SpecialIncoherentResonanceCalculator(int pIdx, unsigned int res_idx);
    EXEC_TARGET devcomplex<fptype> operator()(thrust::tuple<int, fptype*, int> t) const;

private:

    unsigned int resonance_i;
    unsigned int parameters;
};





#endif

