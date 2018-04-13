#pragma once

#include <goofit/PDFs/GooPdf.h>
#include <goofit/detail/Complex.h>

namespace GooFit {

typedef fpcomplex (*resonance_function_ptr)(fptype, fptype, fptype, ParameterContainer &pc);

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

/// Service class intended to hold parametrisations of
/// resonances on Dalitz plots. Don't try to use this
/// as a standalone PDF! It should only be used as a
/// component in one of the friend classes. It extends
/// GooPdf so as to take advantage of the
/// infrastructure, but will crash if used on its own.
class ResonancePdf : public GooPdf {
    friend class TddpPdf;
    friend class DalitzPlotPdf;
    friend class IncoherentSumPdf;

  public:
    ~ResonancePdf() override = default;

    __host__ virtual void recalculateCache() const {}

    void recursiveSetIndices() override;

    __host__ Variable get_amp_real() const { return amp_real; }
    __host__ Variable get_amp_img() const { return amp_imag; }

  protected:
    /// Special constructor that subclasses use
    ResonancePdf(std::string name, Variable ar, Variable ai)
        : GooPdf(name)
        , amp_real(ar)
        , amp_imag(ai) {}

    Variable amp_real;
    Variable amp_imag;

    std::vector<unsigned int> pindices;

    std::vector<fptype> host_constants;

    int resonanceType;
};

namespace Resonances {
/// Relativistic Breit-Wigner
class RBW : public ResonancePdf {
  public:
    RBW(std::string name, Variable ar, Variable ai, Variable mass, Variable width, unsigned int sp, unsigned int cyc);
    ~RBW() override = default;
};

/// LASS
class LASS : public ResonancePdf {
  public:
    LASS(std::string name, Variable ar, Variable ai, Variable mass, Variable width, unsigned int sp, unsigned int cyc);
    ~LASS() override = default;
};

/// Gounaris-Sakurai
class GS : public ResonancePdf {
  public:
    GS(std::string name, Variable ar, Variable ai, Variable mass, Variable width, unsigned int sp, unsigned int cyc);
    ~GS() override = default;
};

/// FLATTE constructor
class FLATTE : public ResonancePdf {
  public:
    FLATTE(std::string name,
           Variable ar,
           Variable ai,
           Variable mean,
           Variable g1,
           Variable rg2og1,
           unsigned int cyc,
           bool symmDP);
    ~FLATTE() override = default;
};

/// Gaussian constructor
class Gauss : public ResonancePdf {
  public:
    Gauss(std::string name, Variable ar, Variable ai, Variable mean, Variable sigma, unsigned int cyc);
    ~Gauss() override = default;
};

/// Nonresonant constructor
class NonRes : public ResonancePdf {
  public:
    NonRes(std::string name, Variable ar, Variable ai);
    ~NonRes() override = default;
};

/// Cubic spline constructor
class Spline : public ResonancePdf {
  public:
    Spline(std::string name,
           Variable ar,
           Variable ai,
           std::vector<fptype> &HH_bin_limits,
           std::vector<Variable> &pwa_coefs_reals,
           std::vector<Variable> &pwa_coefs_imags,
           unsigned int cyc,
           bool symmDP = false);
    ~Spline() override = default;

    /// Recacluate the CACHE values before running
    __host__ void recalculateCache() const override;
};
} // namespace Resonances

} // namespace GooFit
