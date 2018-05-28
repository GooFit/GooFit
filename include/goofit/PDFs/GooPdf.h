#pragma once

#include <thrust/functional.h>

#include <goofit/FitControl.h>
#include <goofit/Log.h>
#include <goofit/PDFs/MetricTaker.h>
#include <goofit/PdfBase.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/Version.h>

#ifdef ROOT_FOUND
class TH1D;
#endif

namespace GooFit {

void *getMetricPointer(EvalFunc val);
//#ifdef SEPARABLE

extern __device__ fptype d_parameters[GOOFIT_MAXPAR];
extern __device__ fptype d_constants[GOOFIT_MAXPAR];
extern __device__ fptype d_observables[GOOFIT_MAXPAR];
extern __device__ fptype d_normalisations[GOOFIT_MAXPAR];

// a couple constants
extern __constant__ fptype c_motherMass;
extern __constant__ fptype c_daug1Mass;
extern __constant__ fptype c_daug2Mass;
extern __constant__ fptype c_daug3Mass;
extern __constant__ fptype c_meson_radius;

// Holds device-side fit parameters.
// extern __constant__ fptype cudaArray[];

// Holds functor-specific indices into cudaArray. Also overloaded to hold integer constants (ie parameters that cannot
// vary.)
// extern __constant__ unsigned int paramIndices[];

// Holds non-integer constants. Notice that first entry is number of events.
// extern __constant__ fptype functorConstants[];

// extern __constant__ fptype normalisationFactors[];

extern __device__ void *device_function_table[GOOFIT_MAXFUNC];
extern void *host_function_table[GOOFIT_MAXFUNC];
extern unsigned int num_device_functions;
extern std::map<void *, int> functionAddressToDeviceIndexMap;
//#endif

// Forward declare
class ParameterContainer;

__device__ int dev_powi(int base, int exp); // Implemented in SmoothHistogramPdf.

/// Pass event, parameters, index into parameters.
typedef fptype (*device_metric_ptr)(fptype, fptype *, fptype);

extern void *host_fcn_ptr;

__device__ fptype callFunction(fptype *eventAddress, ParameterContainer &pc);

class GooPdf : public PdfBase {
  public:
    using PdfBase::PdfBase;

    double calculateNLL() const override;

    /// NB: This does not project correctly in multidimensional datasets, because all observables
    /// other than 'var' will have, for every event, whatever value they happened to get set to last
    /// time they were set. This is likely to be the value from the last event in whatever dataset
    /// you were fitting to, but at any rate you don't get the probability-weighted integral over
    /// the other observables.
    __host__ std::vector<fptype> evaluateAtPoints(Observable var);

    /// A normalize function. This fills in the host_normalize
    __host__ fptype normalize() const override;

    /// Just in case you are British and the previous spelling is offensive
    __host__ fptype normalise() const { return normalize(); }

    __host__ virtual fptype integrate(fptype lo, fptype hi) const { return 0; }
    __host__ bool hasAnalyticIntegral() const override { return false; }
    __host__ fptype getValue(EvalFunc evalfunc = EvalFunc::Eval);

    /// Produce a list of probabilies at points
    __host__ std::vector<std::vector<fptype>> getCompProbsAtDataPoints();

    /// Set an equidistant grid based on the stored variable binning
    __host__ UnbinnedDataSet makeGrid();

    __host__ void initialize();
    __host__ void scan(Observable var, std::vector<fptype> &values);
    __host__ void setFitControl(std::shared_ptr<FitControl> fc) override;
    __host__ virtual void setMetrics();
    __host__ void setParameterConstantness(bool constant = true);

    __host__ virtual void transformGrid(fptype *host_output);
    static __host__ int findFunctionIdx(void *dev_functionPtr);
    __host__ void setDebugMask(int mask, bool setSpecific = true) const;

#ifdef ROOT_FOUND
    /// Plot a PDF to a ROOT histogram
    __host__ TH1D *plotToROOT(Observable var, double normFactor = 1, std::string name = "");
#endif

    __host__ void setIndices() override;

  protected:
    __host__ virtual double sumOfNll(int numVars) const;
    std::shared_ptr<MetricTaker> logger;

  private:
};

} // namespace GooFit
