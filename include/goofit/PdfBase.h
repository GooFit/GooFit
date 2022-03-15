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

// This class never exists on the GPU

class PdfBase {
    friend auto operator<<(std::ostream &, const PdfBase &) -> std::ostream &;

  protected:
    /// Runs once at the beginning of a run. Will always be called, so useful for setup. Inside things like fits, this
    /// will not rerun, however, even if parameters change.
    void pre_run();

    /// This will run before each evaluation (after pre_run), always (even inside fits)
    void pre_call();

    /// use this function to populate the arrays generically, or specialize as needed
    virtual void populateArrays();

    /// This needs to be set before a call to setData.
    void setNumPerTask(PdfBase *p, const int &c);

    /// Generate a range to integrate over
    void generateNormRange();

    // Registration

    /// This adds a parameter.
    void registerParameter(Variable var);

    /// Remove a parameter
    void unregisterParameter(Variable var);

    /// Register a constant
    void registerConstant(fptype value);

    /// Register a function for this PDF to use in evaluation
    template <typename T>
    void registerFunction(std::string name, const T &function) {
        reflex_name_  = name;
        function_ptr_ = get_device_symbol_address(function);
        // add the following 211222 to keep track of the correspondence
        // between function pointer values and their names to make
        // debugging/tracking the executionn of the code easier.  mds
        functionPtrToNameMap[function_ptr_] = name;
    }

    /// Register an observable (Usually done through constructor)
    void registerObservable(Observable obs);

    /// Force all normalization values to 1
    void recursiveSetNormalization(fptype norm = 1.0, bool subpdf = false);

    /// Report the status of the components, observablesList, etc.
    /// added 2101023 mds for debugging
    void status();
    void status(std::string caller);
    // for any Pdf object, print out its name, the name of its class,
    // and the value of the corresponding device function pointer.
    //  added 211026 for debugging
    void listComponents();

  public:
    template <typename... Args>
    explicit PdfBase(std::string pdf_name, std::string n, Args... args)
        : pdf_name_(std::move(pdf_name))
        , name(std::move(n)) {
        std::vector<Observable> obs;
        std::vector<Variable> vars;

        filter_arguments(obs, vars, args...);

        for(auto &ob : obs)
            registerObservable(ob);
        for(auto &var : vars)
            registerParameter(var);
    }

    virtual ~PdfBase() = default;

    // Standard entry functions

    virtual auto calculateNLL() -> double = 0;
    virtual auto normalize() -> fptype    = 0;

    // TODO: Combine with pre_run and remove
    void initializeIndices();

    // TODO: Combine with pre_call and remove
    void copyParams();

    auto getName() const -> std::string { return name; }

    auto getObservables() const -> std::vector<Observable>;
    auto getParameters() const -> std::vector<Variable>;
    auto getParameterByName(std::string n) -> Variable *;

    // User level setup

    void setData(DataSet *data);
    auto getData() -> DataSet * { return data_; }

    virtual void setFitControl(std::shared_ptr<FitControl>) = 0;
    //  for debugging, add a version with an argument that allows
    //  us to track who called this method  mds 211220
    virtual void setFitControl_A(std::shared_ptr<FitControl>, std::string caller) = 0;

    /// Override to indicate that this has an analytic integral
    virtual auto hasAnalyticIntegral() const -> bool { return false; }

    /// Currently only 1D filling supported
    void fillMCDataSimple(size_t events, unsigned int seed = 0);

    /// RooFit style fitting shortcut
    auto fitTo(DataSet *data, int verbosity = 3) -> ROOT::Minuit2::FunctionMinimum;

    /// Even shorter fitting shortcut
    auto fit(int verbosity = 3) -> ROOT::Minuit2::FunctionMinimum;

    auto getFunctionIndex() const -> unsigned int { return functionIdx; }
    auto getParameterIndex() const -> unsigned int { return parameters; }

    void setNormalization(const fptype &v) {
        cachedNormalization                = v;
        host_normalizations[normalIdx + 1] = v;
    }
    auto getNormalization() const -> fptype { return cachedNormalization; }

    /// Set a specific fineness for the integrator
    void setIntegrationFineness(int i);

    /// Have the parameters changed since last evaluation?
    auto parametersChanged() const -> bool;

    void updateVariable(Variable v, fptype newValue);
    void updateParameters();

    // Setup
    void SigGenSetIndices() {
        setupObservables();
        setIndices();
    }
    void setupObservables();
    virtual void recursiveSetIndices();
    virtual void setIndices();

    void setCommonNorm(bool v = true) { commonNorm = v; };
    void setSeparateNorm(bool v = true) { separateNorm = v; };
    auto getCommonNorm() const -> bool { return commonNorm; };
    auto getSeparateNorm() const -> bool { return separateNorm; };
    auto getComponents() -> std::vector<PdfBase *> { return components; };

    /// Get the current PDF name
    auto getPdfName() const -> std::string { return pdf_name_; }

    bool areParamsandConstantsEqualByVal(const PdfBase &other) const;

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

    int integrationBins{-1};        //< Force a specific number of integration bins for all variables
    bool properlyInitialised{true}; //< Allows checking for required extra steps in, eg, Tddp and Convolution.

    unsigned int functionIdx{0}; //< Stores index of device function pointer.

    unsigned int parametersIdx{0};
    unsigned int constantsIdx{0};
    unsigned int observablesIdx{0};
    unsigned int normalIdx{0};

    bool commonNorm{false};
    bool separateNorm{false};

    int m_iEventsPerTask{0};

  private:
    std::string pdf_name_;
    std::string name;
};

auto operator<<(std::ostream &, const PdfBase &) -> std::ostream &;

} // namespace GooFit
