#pragma once

#include <goofit/GlobalCudaDefines.h>

#include <goofit/PDFs/detail/Globals.h>
#include <goofit/Variable.h>
#include <goofit/detail/Abort.h>

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
} // namespace Minuit2
} // namespace ROOT

namespace GooFit {

struct ParameterContainer;

class FitControl;

class DataSet;
class BinnedDataSet;
class UnbinnedDataSet;

namespace {
/// Utility to filter and pick out observables and variables
GOOFIT_MAYBE_UNUSED void filter_arguments(std::vector<Observable> &oblist, std::vector<Variable> &varlist) {}

template <typename... Args>
void filter_arguments(std::vector<Observable> &oblist,
                      std::vector<Variable> &varlist,
                      const Observable &obs,
                      Args... args);
template <typename... Args>
void filter_arguments(std::vector<Observable> &oblist,
                      std::vector<Variable> &varlist,
                      const EventNumber &obs,
                      Args... args);
template <typename... Args>
void filter_arguments(std::vector<Observable> &oblist,
                      std::vector<Variable> &varlist,
                      const Variable &var,
                      Args... args);

/// Pick an Observable off the front
template <typename... Args>
void filter_arguments(std::vector<Observable> &oblist,
                      std::vector<Variable> &varlist,
                      const Observable &obs,
                      Args... args) {
    oblist.push_back(obs);
    return filter_arguments(oblist, varlist, args...);
}

/// Pick an EventNumber off the front
template <typename... Args>
void filter_arguments(std::vector<Observable> &oblist,
                      std::vector<Variable> &varlist,
                      const EventNumber &obs,
                      Args... args) {
    oblist.push_back(obs);
    return filter_arguments(oblist, varlist, args...);
}

/// Pick an Variable (parameter) off the front
template <typename... Args>
void filter_arguments(std::vector<Observable> &oblist,
                      std::vector<Variable> &varlist,
                      const Variable &var,
                      Args... args) {
    varlist.push_back(var);
    return filter_arguments(oblist, varlist, args...);
}
} // namespace

class PdfBase {
    friend std::ostream &operator<<(std::ostream &, const PdfBase &);

  public:
    enum Specials { ForceSeparateNorm = 1, ForceCommonNorm = 2 };

  protected:
    /// Runs once at the beginning of a run. Will always be called, so useful for setup. Inside things like fits, this
    /// will not rerun, however, even if parameters change.
    __host__ void pre_run();

    /// This will run before each evaluation (after pre_run), always (even inside fits)
    __host__ void pre_call();

    /// use this function to populate the arrays generically, or specialize as needed
    virtual void populateArrays();

    /// This needs to be set before a call to setData.
    void setNumPerTask(PdfBase *p, const int &c);

  public:
    template <typename... Args>
    explicit PdfBase(std::string n, Args... args)
        : name(std::move(n)) {
        std::vector<Observable> obs;
        std::vector<Variable> vars;

        filter_arguments(obs, vars, args...);

        for(auto &ob : obs)
            registerObservable(ob);
        for(auto &var : vars)
            registerParameter(var);
    }

    virtual ~PdfBase() = default;

    __host__ virtual double calculateNLL() = 0;
    __host__ virtual fptype normalize()    = 0;
    __host__ void initializeIndices();

    __host__ void addSpecialMask(int m) { specialMask |= m; }
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

    __host__ void fillMCDataSimple(size_t events, unsigned int seed = 0);

    /// RooFit style fitting shortcut
    __host__ ROOT::Minuit2::FunctionMinimum fitTo(DataSet *data, int verbosity = 3);

    /// Even shorter fitting shortcut
    __host__ ROOT::Minuit2::FunctionMinimum fit(int verbosity = 3);

    __host__ unsigned int getFunctionIndex() const { return functionIdx; }
    __host__ unsigned int getParameterIndex() const { return parameters; }

    /// This adds a parameter. The number returned should only be used for checking
    __host__ void registerParameter(Variable var);

    /// The int value returned here is the constant number, for checking
    __host__ void registerConstant(fptype value);

    /// Register a function for this PDF to use in evalution
    template <typename T>
    __host__ void registerFunction(std::string name, const T &function) {
        reflex_name_  = name;
        function_ptr_ = get_device_symbol_address(function);
    }

    /// Register an observable (Usually done through constructor)
    __host__ void registerObservable(Observable obs);

    __host__ void recursiveSetNormalization(fptype norm = 1.0, bool subpdf = false);
    __host__ void unregisterParameter(Variable var);
    __host__ void setIntegrationFineness(int i);

    __host__ bool parametersChanged() const;

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

    std::string reflex_name_; //< This is the name of the type of the PDF, for reflexion purposes. Must be set or
                              // RecursiveSetIndicies must be overloaded.

    void *function_ptr_{nullptr}; //< This is the function pointer to set on the device. Must be set or
                                  // RecursiveSetIndicies must be overloaded.

    fptype numEvents{0};        //< Non-integer to allow weighted events
    unsigned int numEntries{0}; //< Eg number of bins - not always the same as number of events, although it can be.

    fptype *normRanges{nullptr}; //< This is specific to functor instead of variable so that
                                 // MetricTaker::operator needn't use indices.

    unsigned int parameters{0}; //< Stores index, in 'paramIndices', where this functor's information begins.

    fptype cachedNormalization{1.0}; //< Store the normalization for this PDF directly

    std::shared_ptr<FitControl> fitControl;
    std::vector<Observable> observablesList;
    std::vector<Variable> parametersList;
    std::vector<fptype> constantsList;
    std::vector<PdfBase *> components;

    int integrationBins{-1}; //< Force a specific number of integration bins for all variables
    int specialMask{0};      //< For storing information unique to PDFs, eg "Normalize me separately" for TddpPdf.
    bool properlyInitialised{true}; //< Allows checking for required extra steps in, eg, Tddp and Convolution.

    unsigned int functionIdx{0}; //< Stores index of device function pointer.

    unsigned int parametersIdx{0};
    unsigned int constantsIdx{0};
    unsigned int observablesIdx{0};
    unsigned int normalIdx{0};

    int m_iEventsPerTask{0};

  private:
    std::string name;
};

std::ostream &operator<<(std::ostream &, const PdfBase &);

} // namespace GooFit
