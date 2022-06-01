#include <goofit/fitting/FitManagerMinuit2.h>
#include <goofit/fitting/Params.h>

#include <goofit/Color.h>
#include <goofit/Version.h>

#include <goofit/PdfBase.h>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnUserParameters.h>

#include <Minuit2/MnScan.h>
#include <Minuit2/MnMinos.h>

#include <CLI/Timer.hpp>

#include <random>

#ifdef MATHCORE_STANDALONE
#define ROOT_VERSION(x, y, z) 0
#else
#include <RVersion.h>
#endif
#include <iostream>
#include <vector>

namespace GooFit {

FitManagerMinuit2::FitManagerMinuit2(PdfBase *dat)
    : upar_(*dat)
    , pdfPointer(dat)
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

    Minuit2::MnMigrad migrad{fcn_, upar_, strategy};

    // Do the minimization
    if(verbosity > 0)
        std::cout << GooFit::gray << GooFit::bold;

    CLI::Timer avetimer{"Average time per call"};
    Minuit2::FunctionMinimum min = migrad(maxfcn_);

    // Cov Matrix
    matCov = migrad.Covariance();

    if(minos) {
        Minuit2::MnMinos minos{fcn_, min}; // Create MINOS errors
        std::vector<Variable> variables = upar_.GetGooFitParams();
        for(Variable &var : variables) {
            if(var.IsFixed())
                continue;
            else {
                minos_errors.push_back(minos(var.getFitterIndex()));
            }
        }
        // output
        int counter = 0;
        std::cout << "1-sigma minos errors: " << std::endl;
        for(Variable &var : variables) {
            if(var.IsFixed())
                continue;
            else {
                std::cout << var.getName() << ": " << minos_errors[counter].first << " " << minos_errors[counter].second
                          << std::endl;
                counter += 1;
            }
        }
    }

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

auto FitManagerMinuit2::scan(unsigned int par_i, unsigned int maxsteps, fptype low, fptype high)
    -> std::vector<std::pair<fptype, fptype>> {
    std::cout << "\t Scanning Parameter " << upar_.GetName(par_i) << "\n";

    Minuit2::MnScan scan(fcn_, upar_);

    auto vec = scan.Scan(par_i, maxsteps, low, high);

    Minuit2::MnPlot plot;

    plot(vec);

    return vec;
}

void FitManagerMinuit2::printCovMat() { std::cout << std::endl << matCov << std::endl; }

auto FitManagerMinuit2::dmda(fptype a, fptype b) -> fptype {
    fptype ret = a / sqrt(a * a + b * b);
    return ret;
}

auto FitManagerMinuit2::dmdb(fptype a, fptype b) -> fptype {
    fptype ret = b / sqrt(a * a + b * b);
    return ret;
}

auto FitManagerMinuit2::dpda(fptype a, fptype b) -> fptype {
    fptype ret = (-b / (a * a + b * b));
    return ret;
}

auto FitManagerMinuit2::dpdb(fptype a, fptype b) -> fptype {
    fptype ret = (a / (a * a + b * b));
    return ret;
}

void FitManagerMinuit2::printOriginalParams() {
    auto pdfPointer = getPdf();
    auto vec_vars   = pdfPointer->getParameters();
    std::vector<fptype> floatVarVal, floatVarErr;
    for(auto var : vec_vars) {
        if(var.IsFixed())
            continue;
        std::cout << var.getName() << "\t" << var.getValue() << std::endl;
    }
}

std::vector<std::vector<fptype>> FitManagerMinuit2::printParams() {
    auto pdfPointer                = getPdf();
    std::vector<Variable> vec_vars = pdfPointer->getParameters();
    std::vector<fptype> floatVarVal;
    floatVarVal.clear();

    for(Variable &var : vec_vars) {
        if(var.IsFixed())
            continue;
        floatVarVal.push_back(var.getValue());
    }

    std::vector<fptype> vec_mag, vec_mag_err;
    vec_mag.clear();
    vec_mag_err.clear();
    std::vector<fptype> vec_phi, vec_phi_err;
    vec_phi.clear();
    vec_phi_err.clear();
    std::cout << "Fit Parameters Mag and Phase(Degrees) \n";
    std::cout << "*Note: this fuction runs over all free parameters, if you have other parameters than real and "
                 "imaginary coefs fix them before call this function!  \n";
    std::cout << "free parameter resonance: " << floatVarVal.size() / 2 << std::endl;
    std::cout << std::fixed << std::setprecision(8);
    std::cout << std::fixed << "Resonance\tMagnitude\tError\t\tPhase\t\tError" << std::endl;

    for(int i = 0; i < floatVarVal.size(); i += 2) {
        fptype a = floatVarVal[i];
        fptype b = floatVarVal[i + 1];

        fptype mag = sqrt(a * a + b * b);
        fptype phi = atan(b / a) * 180. / M_PI;

        fptype mag_err = dmda(a, b) * dmda(a, b) * matCov(i, i) + dmdb(a, b) * dmdb(a, b) * matCov(i + 1, i + 1)
                         + 2 * dmda(a, b) * dmdb(a, b) * matCov(i, i + 1);

        if(mag_err < 0)
            mag_err = 0;

        mag_err = sqrt(mag_err);

        fptype phi_err = dpda(a, b) * dpda(a, b) * matCov(i, i) + dpdb(a, b) * dpdb(a, b) * matCov(i + 1, i + 1)
                         + 2 * dpda(a, b) * dpdb(a, b) * matCov(i, i + 1);
        if(phi_err < 0)
            phi_err = 0;
        phi_err = sqrt(phi_err) * 180. / M_PI;

        if(a < 0 && b < 0)
            phi -= 180;
        if(a < 0 && b > 0)
            phi += 180;
        vec_mag.push_back(mag);
        vec_phi.push_back(phi);
        vec_mag_err.push_back(mag_err);
        vec_phi_err.push_back(phi_err);
        std::cout << std::fixed << "Resonance " << (i + 2) / 2 << "\t" << mag << "\t" << mag_err << "\t" << phi << "\t"
                  << phi_err << std::endl;
    }

    std::cout << "" << std::endl;

    std::vector<std::vector<fptype>> ret;
    ret.clear();
    ret.push_back(vec_phi);
    ret.push_back(vec_phi_err);
    return ret;
}

void FitManagerMinuit2::setRandMinuitValues(size_t nSamples) {
    std::random_device rd{};
    std::mt19937 gen{rd()};

    std::vector<fptype> floatVarVal;
    floatVarVal.clear();
    std::vector<fptype> floatVarErr;
    floatVarErr.clear();
    auto pdfPointer                = getPdf();
    std::vector<Variable> vec_vars = pdfPointer->getParameters();
    for(Variable &var : vec_vars) {
        if(var.IsFixed())
            continue;
        floatVarVal.push_back(var.getValue());
        floatVarErr.push_back(var.getError());
    }
    const int nFPars = floatVarVal.size();

    VectorXd vy(nFPars);
    samples.clear();

    for(int ii = 0; ii < nSamples; ii++) {
        for(int i = 0; i < nFPars; i++) {
            std::normal_distribution<> d{floatVarVal[i], 1};
            vy(i) = d(gen);
        }

        samples.emplace_back(vy);
    }
}

void FitManagerMinuit2::loadSample(size_t iSample) {
    auto pdfPointer = getPdf();
    auto var        = pdfPointer->getParameters();
    int counter     = 0;
    for(int i = 0; i < var.size(); ++i) {
        if(var[i].IsFixed())
            continue;
        pdfPointer->updateVariable(var[i], samples[iSample](counter));
        counter++;
    }

    Minuit2::MnScan FitManagerMinuit2::getMnScan() {
        Minuit2::MnScan mnscan{fcn_, upar_};

        return mnscan;
    }

} // namespace GooFit
