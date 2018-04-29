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

const int maxParams = GOOFIT_MAXPAR;
extern fptype *dev_event_array;

extern fptype host_parameters[maxParams];
extern fptype host_constants[maxParams];
extern fptype host_observables[maxParams];
extern fptype host_normalisations[maxParams];

extern __constant__ unsigned int c_totalEvents;
extern int totalParameters;
extern int totalConstants;
extern int totalObservables;
extern int totalNormalisations;

class FitControl;

class DataSet;
class BinnedDataSet;
class UnbinnedDataSet;

namespace {
/// Utility to filter and pick out observables and variables
[[ maybe_unused, gnu::unused ]] void filter_arguments(std::vector<Observable> &oblist) {}

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
    __host__ void initializeIndices();

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
    __host__ DataSet *getData() { return data_; }

    __host__ virtual void setFitControl(std::shared_ptr<FitControl>) = 0;
    __host__ virtual bool hasAnalyticIntegral() const { return false; }

    /// RooFit style fitting shortcut
    __host__ ROOT::Minuit2::FunctionMinimum fitTo(DataSet *data, int verbosity = 3);

    __host__ unsigned int getFunctionIndex() const { return functionIdx; }
    __host__ unsigned int getParameterIndex() const { return parameters; }
    // TODO: Make the return value void to find and cleanup old style access
    /// This adds a parameter. The number returned should only be used for checking
    __host__ unsigned int registerParameter(Variable var);
    /// The int value returned here is the constant number, for checking
    __host__ unsigned int registerConstant(fptype value);
    __host__ unsigned int registerConstants(unsigned int amount);
    __host__ virtual void recursiveSetNormalisation(fptype norm = 1) const;
    __host__ void unregisterParameter(Variable var);
    __host__ void registerObservable(Observable obs);
    __host__ void setIntegrationFineness(int i);
    __host__ void printProfileInfo(bool topLevel = true);

    __host__ bool parametersChanged() const;

    __host__ void checkInitStatus(std::vector<std::string> &unInited) const;
    __host__ void clearCurrentFit();
    __host__ void SigGenSetIndices() {
        setupObservables();
        setIndices();
    }

    __host__ void updateVariable(Variable v, fptype newValue);
    __host__ void updateParameters();

    __host__ void setupObservables();

    __host__ virtual void recursiveSetIndices();

    __host__ virtual void setIndices();

  protected:
    DataSet *data_ = nullptr; //< Remember the original dataset
    // use this function to populate the arrays generically.
    __host__ void populateArrays();

    fptype numEvents{0};        //< Non-integer to allow weighted events
    unsigned int numEntries{0}; //< Eg number of bins - not always the same as number of events, although it can be.
    fptype *normRanges{
        nullptr}; //< This is specific to functor instead of variable so that MetricTaker::operator needn't use indices.
    unsigned int parameters{0}; //< Stores index, in 'paramIndices', where this functor's information begins.
    unsigned int cIndex{1};     //< Stores location of constants.
    std::vector<Observable> observablesList;
    std::vector<Variable> parametersList;
    std::vector<fptype> constantsList;
    std::shared_ptr<FitControl> fitControl;
    std::vector<PdfBase *> components;
    int integrationBins{-1};
    int specialMask{0}; //< For storing information unique to PDFs, eg "Normalize me separately" for TddpPdf.
    bool properlyInitialised{true}; //< Allows checking for required extra steps in, eg, Tddp and Convolution.

    unsigned int functionIdx{0}; //< Stores index of device function pointer.

    unsigned int parametersIdx{0};
    unsigned int constantsIdx{0};
    unsigned int observablesIdx{0};
    unsigned int normalIdx{0};
    int m_iEventsPerTask{0};

    /// This needs to be set before a call to setData.
    void setNumPerTask(PdfBase *p, const int &c);

  private:
    std::string name;
};

} // namespace GooFit
