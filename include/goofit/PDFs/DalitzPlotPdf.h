#ifndef DALITZPLOT_PDF_HH
#define DALITZPLOT_PDF_HH

#include "goofit/PDFs/GooPdf.h"
#include "goofit/PDFs/DalitzPlotHelpers.h"

#include <thrust/complex.h>

class SpecialResonanceIntegrator;
class SpecialResonanceCalculator;

class DalitzPlotPdf : public GooPdf {
public:
    DalitzPlotPdf(std::string n, Variable* m12, Variable* m13, CountingVariable* eventNumber, DecayInfo* decay,
                  GooPdf* eff);
    // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the
    // coherent sum. The caching method requires that it be done this way or the ProdPdf
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
    DEVICE_VECTOR<thrust::complex<fptype>>* cachedWaves[16]; // Caches the BW values for each event.
    thrust::complex<fptype>*** integrals; // Caches the integrals of the BW waves for each combination of resonances.

    bool* redoIntegral;
    mutable bool forceRedoIntegrals;
    fptype* cachedMasses;
    fptype* cachedWidths;
    int totalEventSize;
    int cacheToUse;
    SpecialResonanceIntegrator*** integrators;
    SpecialResonanceCalculator** calculators;
};

class SpecialResonanceIntegrator : public thrust::unary_function<thrust::tuple<int, fptype*>, thrust::complex<fptype>> {
public:
    // Class used to calculate integrals of terms BW_i * BW_j^*.
    SpecialResonanceIntegrator(int pIdx, unsigned int ri, unsigned int rj);
    __device__ thrust::complex<fptype> operator()(thrust::tuple<int, fptype*> t) const;
private:

    unsigned int resonance_i;
    unsigned int resonance_j;
    unsigned int parameters;
};

class SpecialResonanceCalculator : public thrust::unary_function<thrust::tuple<int, fptype*, int>, thrust::complex<fptype>> {
public:
    // Used to create the cached BW values.
    SpecialResonanceCalculator(int pIdx, unsigned int res_idx);
    __device__ thrust::complex<fptype> operator()(thrust::tuple<int, fptype*, int> t) const;

private:

    unsigned int resonance_i;
    unsigned int parameters;
};


#endif

