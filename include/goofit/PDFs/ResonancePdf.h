#ifndef RESONANCE_PDF_HH
#define RESONANCE_PDF_HH

#include "goofit/PDFs/GooPdf.h"
#include "goofit/PDFs/devcomplex.h"
typedef devcomplex<fptype> (*resonance_function_ptr)(fptype, fptype, fptype, unsigned int*);

__device__ fptype twoBodyCMmom(double rMassSq, fptype d1m, fptype d2m);

__device__ fptype dampingFactorSquare(const fptype& cmmom, const int& spin, const fptype& mRadius);

__device__ fptype spinFactor(unsigned int spin, fptype motherMass, fptype daug1Mass, fptype daug2Mass,
                              fptype daug3Mass, fptype m12, fptype m13, fptype m23, unsigned int cyclic_index);

class ResonancePdf : public GooPdf {
    // Service class intended to hold parametrisations of
    // resonances on Dalitz plots. Don't try to use this
    // as a standalone PDF! It should only be used as a
    // component in one of the friend classes. It extends
    // GooPdf so as to take advantage of the
    // infrastructure, but will crash if used on its own.

    friend class TddpPdf;
    friend class DalitzPlotPdf;
    friend class IncoherentSumPdf;
public:
    // Constructor for regular BW
    ResonancePdf(string name,
                 Variable* ar,
                 Variable* ai,
                 Variable* mass,
                 Variable* width,
                 unsigned int sp,
                 unsigned int cyc);

    // Gounaris-Sakurai
    ResonancePdf(string name,
                 Variable* ar,
                 Variable* ai,
                 unsigned int sp,
                 Variable* mass,
                 Variable* width,
                 unsigned int cyc);

    // LASS constructor
    ResonancePdf(string name,
                 Variable* ar,
                 Variable* ai,
                 Variable* mass,
                 unsigned int sp,
                 Variable* width,
                 unsigned int cyc);


    // Nonresonant constructor
    ResonancePdf(string name,
                 Variable* ar,
                 Variable* ai);

    // Gaussian constructor
    ResonancePdf(string name,
                 Variable* ar,
                 Variable* ai,
                 Variable* mean,
                 Variable* sigma,
                 unsigned int cyc);

private:
    void setConstantIndex(unsigned int idx) {
        host_indices[parameters + 1] = idx;
    }

    Variable* amp_real;
    Variable* amp_imag;
    /*
    Variable* mass;
    Variable* width;
    unsigned int spin;
    unsigned int cyclic_index;
    unsigned int eval_type;
    unsigned int resonance_type;
    */
};

#endif
