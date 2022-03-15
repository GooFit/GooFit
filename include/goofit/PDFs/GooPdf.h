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
__device__ auto callFunction(fptype *eventAddress, ParameterContainer &pc) -> fptype;

class GooPdf : public PdfBase {
  protected:
    // These are helper methods that do all the work

    /// This collects the number of variables for the thrust call
    /// -(n+2) for binned evalute's
    auto get_event_size() const -> int;

    /// This reduces the current function over the data. Does *not* prepare
    /// or normalize TODO: Remove const if caching added
    __host__ auto reduce_with_metric() const -> double;

    /// This evaluates the current function over the data. Does *not* prepare
    /// or normalize  TODO: Remove const if caching added
    __host__ void evaluate_with_metric(thrust::device_vector<fptype> &results) const;

    /// Shortcut to allow .cpp files to use this if they need to
    __host__ auto evaluate_with_metric() const -> thrust::host_vector<fptype>;

  public:
    using PdfBase::PdfBase;

    ~GooPdf() override;

    /// This is the total number of bins in the normalize grid
    __host__ auto get_bin_grid_size() const -> int;

    /// This is the total ND Volume of a cube in the normalize grid
    __host__ auto get_bin_grid_volume() const -> fptype;

    /// This reduces the current function over the bin grid. Does *not* prepare
    /// or normalize. Used by normalize.  TODO: Remove const if caching added
    __host__ auto reduce_with_bins() const -> double;

    __host__ auto calculateNLL() -> double override;

    /// NB: This does not project correctly in multidimensional datasets, because all observables
    /// other than 'var' will have, for every event, whatever value they happened to get set to last
    /// time they were set. This is likely to be the value from the last event in whatever dataset
    /// you were fitting to, but at any rate you don't get the probability-weighted integral over
    /// the other observables.
    __host__ auto evaluateAtPoints(Observable var) -> std::vector<fptype>;

    /// A normalize function. This fills in the host_normalize
    __host__ auto normalize() -> fptype override;

    /// Just in case you are British and the previous spelling is offensive
    __host__ auto normalise() -> fptype { return normalize(); }

    __host__ virtual auto integrate(fptype lo, fptype hi) const -> fptype;
    __host__ auto hasAnalyticIntegral() const -> bool override;
    __host__ auto getValue(EvalFunc evalfunc = EvalFunc::Eval) -> fptype;

    __host__ void listAllComponents(PdfBase *someComponent);

    /// Produce a list of probabilies at points
    __host__ auto getCompProbsAtDataPoints() -> std::vector<std::vector<fptype>>;

    /// Set an equidistant grid based on the stored variable binning
    __host__ auto makeGrid() -> UnbinnedDataSet;

    __host__ void initialize();
    __host__ void scan(Observable var, std::vector<fptype> &values);
    __host__ void setFitControl(std::shared_ptr<FitControl> fc) override;
    //  for debugging, add a version with an argument that allows
    //  us to track who called this method  mds 211220
    __host__ void setFitControl_A(std::shared_ptr<FitControl>, std::string caller) override;
    __host__ virtual void setMetrics();
    __host__ void setParameterConstantness(bool constant = true);

    static __host__ auto findFunctionIdx(void *dev_functionPtr) -> int;
    static __host__ auto lookUpFunctionIdx(void *dev_functionPtr) -> int;
    __host__ void setDebugMask(int mask, bool setSpecific = true) const;

#ifdef ROOT_FOUND
    /// Plot a PDF to a ROOT histogram
    __host__ auto plotToROOT(Observable var, double normFactor = 1, std::string name = "") -> TH1D *;
#endif

    __host__ void setIndices() override;

  protected:
    std::shared_ptr<MetricTaker> logger;
};

} // namespace GooFit
