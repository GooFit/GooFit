#pragma once

#include <goofit/FitControl.h>
#include <goofit/PDFs/MetricPointer.h>
#include <goofit/PDFs/MetricTaker.h>
#include <goofit/PDFs/detail/Globals.h>
#include <goofit/PdfBase.h>
#include <goofit/UnbinnedDataSet.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#ifdef ROOT_FOUND
class TH1D;
#endif

namespace GooFit {

// Forward declare
struct ParameterContainer;
__device__ fptype callFunction(fptype *eventAddress, ParameterContainer &pc);

class GooPdf : public PdfBase {
  protected:
    // These are helper methods that do all the work

    /// This collects the number of variables for the thrust call
    /// -(n+2) for binned evalutes
    int get_event_size() const;

    /// This is the total number of bins in the normalize grid
    int get_bin_grid_size() const;

    /// This is the total ND Volume of a cube in the normalize grid
    fptype get_bin_grid_volume() const;

    /// This reduces the current function over the data. Does *not* prepare
    /// or normalize TODO: Remove const if caching added
    __host__ double reduce_with_metric() const;

    /// This reduces the current function over the bin grid. Does *not* prepare
    /// or normalize. Used by normalize.  TODO: Remove const if caching added
    __host__ double reduce_with_bins() const;

    /// This evaluates the current function over the data. Does *not* prepare
    /// or normalize  TODO: Remove const if caching added
    __host__ void evaluate_with_metric(thrust::device_vector<fptype> &results) const;

    /// Shortcut to allow .cpp files to use this if they need to
    __host__ thrust::host_vector<fptype> evaluate_with_metric() const;

  public:
    using PdfBase::PdfBase;

    ~GooPdf() override;

    __host__ double calculateNLL() override;

    /// NB: This does not project correctly in multidimensional datasets, because all observables
    /// other than 'var' will have, for every event, whatever value they happened to get set to last
    /// time they were set. This is likely to be the value from the last event in whatever dataset
    /// you were fitting to, but at any rate you don't get the probability-weighted integral over
    /// the other observables.
    __host__ std::vector<fptype> evaluateAtPoints(Observable var);

    /// A normalize function. This fills in the host_normalize
    __host__ fptype normalize() override;

    /// Just in case you are British and the previous spelling is offensive
    __host__ fptype normalise() { return normalize(); }

    __host__ virtual fptype integrate(fptype lo, fptype hi) const;
    __host__ bool hasAnalyticIntegral() const override;
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

    static __host__ int findFunctionIdx(void *dev_functionPtr);
    __host__ void setDebugMask(int mask, bool setSpecific = true) const;

#ifdef ROOT_FOUND
    /// Plot a PDF to a ROOT histogram
    __host__ TH1D *plotToROOT(Observable var, double normFactor = 1, std::string name = "");
#endif

    __host__ void setIndices() override;

  protected:
    std::shared_ptr<MetricTaker> logger;
};

} // namespace GooFit
