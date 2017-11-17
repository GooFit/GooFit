#pragma once

#include <goofit/GlobalCudaDefines.h>

#include <goofit/Variable.h>
#include <goofit/Version.h>
#include <goofit/detail/Abort.h>
#include <goofit/detail/Macros.h>

#include <algorithm>
#include <map>
#include <set>
#include <vector>

#ifdef ROOT_FOUND
class TH1D;
#endif

namespace ROOT {
namespace Minuit2 {
class FunctionMinimum;
}
} // namespace ROOT

namespace GooFit {

/* Future use, apperently:
#include <thrust/device_vector.h>
#include <thrust/iterator/constant_iterator.h>

typedef thrust::counting_iterator<int> IndexIterator;
typedef thrust::constant_iterator<fptype*> DataIterator;
typedef thrust::constant_iterator<int> SizeIterator;
typedef thrust::tuple<IndexIterator, DataIterator, SizeIterator> EventTuple;
typedef thrust::zip_iterator<EventTuple> EventIterator;
*/

const int maxParams = GOOFIT_MAXPAR;
extern fptype *dev_event_array;
extern fptype host_normalisation[maxParams];
extern fptype host_params[maxParams];
extern unsigned int host_indices[maxParams];
extern int totalParams;
extern int totalConstants;

class FitControl;

class DataSet;
class BinnedDataSet;
class UnbinnedDataSet;

namespace {
/// Utility to filter and pick out observables and variables
void filter_arguments(std::vector<Observable> &oblist) {}

template <typename... Args>
void filter_arguments(std::vector<Observable> &oblist, const Observable &obs, Args... args) {
    oblist.push_back(obs);
    return filter_arguments(oblist, args...);
}

template <typename... Args>
void filter_arguments(std::vector<Observable> &oblist, const EventNumber &obs, Args... args) {
    oblist.push_back(obs);
    return filter_arguments(oblist, args...);
}
} // namespace

class PdfBase {
  public:
    template <typename... Args>
    explicit PdfBase(std::string n, Args... args)
        : name(std::move(n)) {
        std::vector<Observable> obs;
        filter_arguments(obs, args...);
        for(auto &ob : obs)
            registerObservable(ob);
    }

    virtual ~PdfBase() = default;

    enum Specials { ForceSeparateNorm = 1, ForceCommonNorm = 2 };

    __host__ virtual double calculateNLL() const = 0;
    __host__ virtual fptype normalize() const    = 0;
    __host__ void initializeIndices(std::vector<unsigned int> pindices);

    __host__ void addSpecialMask(int m) { specialMask |= m; }
    __host__ void copyParams(const std::vector<double> &pars) const;
    __host__ void copyParams();
    __host__ void copyNormFactors() const;
    __host__ void generateNormRange();
    __host__ std::string getName() const { return name; }

    __host__ virtual std::vector<Observable> getObservables() const;
    __host__ virtual std::vector<Variable> getParameters() const;

    __host__ Variable *getParameterByName(std::string n);
    __host__ int getSpecialMask() const { return specialMask; }

    __host__ void setData(DataSet *data);
    __host__ DataSet *getData();

    __host__ virtual void setFitControl(std::shared_ptr<FitControl>) = 0;
    __host__ virtual bool hasAnalyticIntegral() const { return false; }

    /// RooFit style fitting shortcut
    __host__ ROOT::Minuit2::FunctionMinimum fitTo(DataSet *data, int verbosity = 3);

    __host__ unsigned int getFunctionIndex() const { return functionIdx; }
    __host__ unsigned int getParameterIndex() const { return parameters; }
    __host__ unsigned int registerParameter(Variable var);
    __host__ unsigned int registerConstants(unsigned int amount);
    __host__ virtual void recursiveSetNormalisation(fptype norm = 1) const;
    __host__ void unregisterParameter(Variable var);
    __host__ void registerObservable(Observable obs);
    __host__ void setIntegrationFineness(int i);
    __host__ void printProfileInfo(bool topLevel = true);

    __host__ bool parametersChanged() const;

    __host__ void checkInitStatus(std::vector<std::string> &unInited) const;
    void clearCurrentFit();
    __host__ void SigGenSetIndices() { setIndices(); }

  protected:
    DataSet *data_ = nullptr;   //< Remember the original dataset
    fptype numEvents{0};        //< Non-integer to allow weighted events
    unsigned int numEntries{0}; //< Eg number of bins - not always the same as number of events, although it can be.
    fptype *normRanges{
        nullptr}; //< This is specific to functor instead of variable so that MetricTaker::operator needn't use indices.
    unsigned int parameters{0}; //< Stores index, in 'paramIndices', where this functor's information begins.
    unsigned int cIndex{1};     //< Stores location of constants.
    std::vector<Observable> observables;
    std::vector<Variable> parameterList;
    std::shared_ptr<FitControl> fitControl;
    std::vector<PdfBase *> components;
    int integrationBins{-1};
    int specialMask{0}; //< For storing information unique to PDFs, eg "Normalize me separately" for TddpPdf.
    bool properlyInitialised{true}; //< Allows checking for required extra steps in, eg, Tddp and Convolution.

    unsigned int functionIdx{0}; //< Stores index of device function pointer.

    int m_iEventsPerTask{0};

    /// This needs to be set before a call to setData.
    void setNumPerTask(PdfBase *p, const int &c);

  private:
    std::string name;

    __host__ void recursiveSetIndices();
    __host__ void setIndices();
};

} // namespace GooFit
