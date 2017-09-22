#pragma once

#include "goofit/PDFs/GooPdf.h"
#include "goofit/detail/Complex.h"

namespace GooFit {

enum class ResPdfType { RBW = 0, LASS, GS, FLATTE, GAUSS, SPLINE, NONRES };

#define MAXNKNOBS 1000

typedef fpcomplex (*resonance_function_ptr)(fptype, fptype, fptype, unsigned int *);

__device__ fptype twoBodyCMmom(double rMassSq, fptype d1m, fptype d2m);

__device__ fptype dampingFactorSquare(const fptype &cmmom, const int &spin, const fptype &mRadius);

__device__ fptype spinFactor(unsigned int spin,
                             fptype motherMass,
                             fptype daug1Mass,
                             fptype daug2Mass,
                             fptype daug3Mass,
                             fptype m12,
                             fptype m13,
                             fptype m23,
                             unsigned int cyclic_index);

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
    // Constructor for regular BW,Gounaris-Sakurai,LASS
    ResonancePdf(std::string name,
                 ResPdfType rpt,
                 Variable *ar,
                 Variable *ai,
                 Variable *mass,
                 Variable *width,
                 unsigned int sp,
                 unsigned int cyc,
                 bool symmDP = false);

    // Constructor for NONRES
    ResonancePdf(std::string name, ResPdfType rpt, Variable *ar, Variable *ai);

    // Gaussian constructor
    ResonancePdf(std::string name,
                 ResPdfType rpt,
                 Variable *ar,
                 Variable *ai,
                 Variable *mean,
                 Variable *sigma,
                 unsigned int cyc);

    // Flatte constructor (arXiv:1505.01710)
    ResonancePdf(std::string name,
                 ResPdfType rpt,
                 Variable *ar,
                 Variable *ai,
                 Variable *mean,
                 Variable *g1,
                 Variable *rg2og1,
                 unsigned int cyc,
                 const bool symmDP);

    // Cubic spline constructor
    ResonancePdf(std::string name,
                 ResPdfType rpt,
                 Variable *ar,
                 Variable *ai,
                 std::vector<fptype> &HH_bin_limits,
                 std::vector<Variable *> &pwa_coefs_reals,
                 std::vector<Variable *> &pwa_coefs_imags,
                 unsigned int cyc,
                 const bool symmDP = false);

    __host__ void recalculateCache() const;

  private:
    void setConstantIndex(unsigned int idx) { host_indices[parameters + 1] = idx; }

    Variable *amp_real;
    Variable *amp_imag;
    /*
    Variable* mass;
    Variable* width;
    unsigned int spin;
    unsigned int cyclic_index;
    unsigned int eval_type;
    unsigned int resonance_type;
    */

    std::vector<fptype> host_constants;

    const ResPdfType rpt_;
};

} // namespace GooFit
