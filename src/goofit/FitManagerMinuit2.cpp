#include <goofit/fitting/FitManagerMinuit2.h>

#include <goofit/Color.h>
#include <goofit/Version.h>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnUserParameters.h>

#include <CLI/Timer.hpp>

#ifdef MATHCORE_STANDALONE
#define ROOT_VERSION(x, y, z) 0
#else
#include <RVersion.h>
#endif

namespace GooFit {

FitManagerMinuit2::FitManagerMinuit2(PdfBase *dat)
    : upar_(*dat)
    , fcn_(upar_) {}

auto FitManagerMinuit2::fit() -> Minuit2::FunctionMinimum {
#if !defined(MATHCORE_STANDALONE) && GOOFIT_ROOT_FOUND && ROOT_VERSION_CODE < ROOT_VERSION(6, 24, 0)
    auto val = Minuit2::MnPrint::Level();
    Minuit2::MnPrint::SetLevel(verbosity);
#else
    auto val = Minuit2::MnPrint::GlobalLevel();
    Minuit2::MnPrint::SetGlobalLevel(verbosity);
#endif

    // Setting global call number to 0
    host_callnumber = 0;

    CLI::Timer timer{"The minimization took"};

    Minuit2::MnMigrad migrad{fcn_, upar_};

    // Do the minimization
    if(verbosity > 0)
        std::cout << GooFit::gray << GooFit::bold;

    CLI::Timer avetimer{"Average time per call"};
    Minuit2::FunctionMinimum min = migrad(maxfcn_);

    // Print nice output
    if(verbosity > 0) {
        std::cout << GooFit::reset << (min.IsValid() ? GooFit::green : GooFit::red);
        std::cout << min << GooFit::reset;
        std::cout << GooFit::magenta << timer.to_string() << std::endl;
        std::cout << (avetimer / min.NFcn()).to_string() << GooFit::reset << std::endl;
    }

    if(min.IsValid() && min.HasCovariance() && !min.IsAboveMaxEdm() && !min.HasReachedCallLimit()) {
        retval_ = FitErrors::Valid;
    } else {
        if(verbosity > 0) {
            std::cout << GooFit::red;
            std::cout << "HesseFailed: " << min.HesseFailed() << std::endl;
            std::cout << "HasCovariance: " << min.HasCovariance() << std::endl;
            std::cout << "HasValidCovariance: " << min.HasValidCovariance() << std::endl;
            std::cout << "HasValidParameters: " << min.HasValidParameters() << std::endl;
            std::cout << "IsAboveMaxEdm: " << min.IsAboveMaxEdm() << std::endl;
            std::cout << "HasReachedCallLimit: " << min.HasReachedCallLimit() << std::endl;
            std::cout << "HasAccurateCovar: " << min.HasAccurateCovar() << std::endl;
            std::cout << "HasPosDefCovar : " << min.HasPosDefCovar() << std::endl;
            std::cout << "HasMadePosDefCovar : " << min.HasMadePosDefCovar() << std::endl;
            std::cout << GooFit::reset;
        }

        retval_ = FitErrors::InValid;
    }

    // Set the parameters in GooFit to the new values
    upar_.SetGooFitParams(min.UserState());

#if !defined(MATHCORE_STANDALONE) && GOOFIT_ROOT_FOUND && ROOT_VERSION_CODE < ROOT_VERSION(6, 24, 0)
    Minuit2::MnPrint::SetLevel(val);
#else
    Minuit2::MnPrint::SetGlobalLevel(val);
#endif
    return min;
}

} // namespace GooFit
