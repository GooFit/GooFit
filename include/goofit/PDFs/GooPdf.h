#pragma once

#include <thrust/functional.h>

#include "goofit/PDFs/MetricTaker.h"
#include "goofit/PdfBase.h"
#include "goofit/UnbinnedDataSet.h"
#include "goofit/PDFs/MetricTaker.h"
#include "goofit/Log.h"

#ifdef ROOT_FOUND
class TH1D;
#endif

namespace GooFit {

// TODO: Replace this with class MetricTaker;
// And fill in the .cu files where needed

#ifdef SEPARABLE

struct ParameterContainer {
    __host__ __device__ ParameterContainer();
    __host__ __device__ ParameterContainer(const ParameterContainer &pc);

    fptype *parameters;
    fptype *constants;
    fptype *observables;
    fptype *normalisations;

    int parameterIdx;
    int constantIdx;
    int observableIdx;
    int normalIdx;

    int funcIdx;

    // each PDF needs to supply the amount of each array used.
    // This function automatically adds +1 for the size.
    __device__ void incrementIndex(const int funcs, const int params, const int cons, const int obs, const int norms) {
        funcIdx += funcs;
        parameterIdx += params + 1;
        constantIdx += cons + 1;
        observableIdx += obs + 1;
        normalIdx += norms + 1;
    }

    // slow version, avoid at all costs!
    __device__ void incrementIndex() {
        funcIdx++;

        int np = parameters[parameterIdx];
        int nc = constants[constantIdx];
        int no = observables[observableIdx];
        int nn = normalisations[normalIdx];

        parameterIdx += np + 1;
        constantIdx += nc + 1;
        observableIdx += no + 1;
        normalIdx += nn + 1;
    }
};

extern __device__ fptype d_parameters[maxParams];
extern __device__ fptype d_constants[maxParams];
extern __device__ fptype d_observables[maxParams];
extern __device__ fptype d_normalisations[maxParams];

// a couple constants
extern __constant__ fptype c_motherMass;
extern __constant__ fptype c_daug1Mass;
extern __constant__ fptype c_daug2Mass;
extern __constant__ fptype c_daug3Mass;
extern __constant__ fptype c_meson_radius;

// Holds device-side fit parameters.
// extern __constant__ fptype cudaArray[maxParams];

// Holds functor-specific indices into cudaArray. Also overloaded to hold integer constants (ie parameters that cannot
// vary.)
// extern __constant__ unsigned int paramIndices[maxParams];

// Holds non-integer constants. Notice that first entry is number of events.
// extern __constant__ fptype functorConstants[maxParams];

// extern __constant__ fptype normalisationFactors[maxParams];

extern __device__ void *device_function_table[200];
extern void *host_function_table[200];
extern unsigned int num_device_functions;
extern std::map<void *, int> functionAddressToDeviceIndexMap;
#endif

__device__ int dev_powi(int base, int exp); // Implemented in SmoothHistogramPdf.
void *getMetricPointer(std::string name);

/// Pass event, parameters, index into parameters.
typedef fptype (*device_function_ptr)(fptype *, ParameterContainer &);

typedef fptype (*device_metric_ptr)(fptype, fptype, fptype);

extern void *host_fcn_ptr;

__device__ fptype callFunction(fptype *eventAddress, ParameterContainer &pc);

class GooPdf : public PdfBase {
  public:
    GooPdf(Variable *x, std::string n);
    __host__ double calculateNLL() const override;

    /// NB: This does not project correctly in multidimensional datasets, because all observables
    /// other than 'var' will have, for every event, whatever value they happened to get set to last
    /// time they were set. This is likely to be the value from the last event in whatever dataset
    /// you were fitting to, but at any rate you don't get the probability-weighted integral over
    /// the other observables.
    __host__ std::vector<fptype> evaluateAtPoints(Variable *var);

    /// A normalize function. This fills in the host_normalize
    __host__ fptype normalize() const override;

    /// Just in case you are British and the previous spelling is offensive
    __host__ fptype normalise() const { return normalize(); }

    __host__ virtual fptype integrate(fptype lo, fptype hi) const { return 0; }
    __host__ bool hasAnalyticIntegral() const override { return false; }
    __host__ fptype getValue();

    /// Produce a list of probabilies at points
    __host__ std::vector<std::vector<fptype>> getCompProbsAtDataPoints();

    /// Set an equidistant grid based on the stored variable binning
    __host__ UnbinnedDataSet makeGrid();

    __host__ void initialize(std::vector<unsigned int> pindices, void *dev_functionPtr = host_fcn_ptr);
    __host__ void scan(Variable *var, std::vector<fptype> &values);
    __host__ void setFitControl(FitControl *const fc, bool takeOwnerShip = true) override;
    __host__ virtual void setMetrics();
    __host__ void setParameterConstantness(bool constant = true);

    __host__ virtual void transformGrid(fptype *host_output);
    static __host__ int findFunctionIdx(void *dev_functionPtr);
    __host__ void setDebugMask(int mask, bool setSpecific = true) const;

#ifdef ROOT_FOUND
    /// Plot a PDF to a ROOT histogram
    __host__ TH1D *plotToROOT(Variable *var, double normFactor = 1, std::string name = "");
#endif

  protected:
    __host__ virtual void setIndices();

    __host__ virtual double sumOfNll(int numVars) const;
    MetricTaker *logger;

  private:
};

} // namespace GooFit
