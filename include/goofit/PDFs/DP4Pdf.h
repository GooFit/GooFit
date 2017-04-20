/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!
See *.cu file for more details
*/

#ifndef D_P4_PDF_HH
#define D_P4_PDF_HH

#include "goofit/PDFs/GooPdf.h"
#include "goofit/PDFs/DalitzPlotHelpers.h"
#include "goofit/PDFs/SpinFactors.h"
#include <mcbooster/GContainers.h>
#include <tuple>
#include <thrust/remove.h>


#ifdef SEPARABLE
extern __constant__ unsigned int AmpIndices[500];
#endif

class LSCalculator;
class AmpCalc;
class SFCalculator;
class NormIntegrator;



class DPPdf : public GooPdf {
public:
    DPPdf(std::string n, std::vector<Variable*> observables, DecayInfo_DP* decay, GooPdf* eff,
          unsigned int MCeventsNorm = 5e6);
    // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the
    // coherent sum. The caching method requires that it be done this way or the ProdPdf
    // normalisation will get *really* confused and give wrong answers.

    __host__ virtual fptype normalise() const;
    __host__ void setDataSize(unsigned int dataSize, unsigned int evtSize = 6);
    __host__ void setForceIntegrals(bool f = true) {
        forceRedoIntegrals = f;
    }
    __host__ int getMCevents() {
        return MCevents;
    }
    __host__ void setGenerationOffset(int off) {
        generation_offset = off;
    }
    __host__ std::tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h,  mcbooster::RealVector_h>
    GenerateSig(unsigned int numEvents);

protected:

private:

    std::map<std::string, std::pair<std::vector<unsigned int>, std::vector<unsigned int>>> AmpMap;
    std::map<std::string, unsigned int> compMap;
    // std::map<unsigned int, unsigned int> massmap;
    std::map<std::string, unsigned int> SpinMap;
    std::vector<SpinFactor*> SpinFactors;
    std::vector<AmpCalc*> AmpCalcs;
    NormIntegrator* Integrator;
    std::vector<SFCalculator*> sfcalculators;
    std::vector<LSCalculator*> lscalculators;

    // store normalization events
    mcbooster::RealVector_d norm_M12;
    mcbooster::RealVector_d norm_M34;
    mcbooster::RealVector_d norm_CosTheta12;
    mcbooster::RealVector_d norm_CosTheta34;
    mcbooster::RealVector_d norm_phi;

    //store spin and lineshape values for normalization
    mutable mcbooster::RealVector_d norm_SF;
    mutable mcbooster::mc_device_vector<devcomplex<fptype>> norm_LS;


    DecayInfo_DP* decayInfo;
    std::vector<Variable*> _observables;
    fptype* hostphsp;
    int MCevents;
    // Following variables are useful if masses and widths, involved in difficult BW calculation,
    // change infrequently while amplitudes, only used in adding BW results together, change rapidly.
    DEVICE_VECTOR<devcomplex<fptype>>* cachedResSF; // Caches the BW values and Spins for each event.
    DEVICE_VECTOR<devcomplex<fptype>>* cachedAMPs; // cache Amplitude values for each event.

    mutable bool generation_no_norm;
    mutable bool SpinsCalculated;
    bool* redoIntegral;
    mutable bool forceRedoIntegrals;
    fptype* cachedMasses;
    fptype* cachedWidths;
    int totalEventSize;
    int cacheToUse;
    int generation_offset;
};


class SFCalculator : public thrust::unary_function<thrust::tuple<int, fptype*, int>, devcomplex<fptype>> {
public:
    // Used to create the cached BW values.
    SFCalculator(int pIdx, unsigned int sf_idx);
    EXEC_TARGET devcomplex<fptype> operator()(thrust::tuple<int, fptype*, int> t) const;

private:

    unsigned int _spinfactor_i;
    unsigned int _parameters;
};

class NormSpinCalculator : public
    thrust::unary_function<thrust::tuple<fptype, fptype, fptype, fptype, fptype>, fptype> {
public:
    // Used to create the cached BW values.
    NormSpinCalculator(int pIdx, unsigned int sf_idx);
    EXEC_TARGET fptype operator()(thrust::tuple<fptype, fptype, fptype, fptype, fptype> t) const;

private:

    unsigned int _spinfactor_i;
    unsigned int _parameters;
};


class LSCalculator : public thrust::unary_function<thrust::tuple<int, fptype*, int>, devcomplex<fptype>> {
public:
    // Used to create the cached BW values.
    LSCalculator(int pIdx, unsigned int res_idx);
    EXEC_TARGET devcomplex<fptype> operator()(thrust::tuple<int, fptype*, int> t) const;

private:

    unsigned int _resonance_i;
    unsigned int _parameters;
};

class NormLSCalculator : public
    thrust::unary_function<thrust::tuple<mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t>, devcomplex<fptype>> {
public:
    // Used to create the cached BW values.
    NormLSCalculator(int pIdx, unsigned int res_idx);
    EXEC_TARGET devcomplex<fptype> operator()(
        thrust::tuple<mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t, mcbooster::GReal_t> t)
    const;

private:

    unsigned int _resonance_i;
    unsigned int _parameters;
};

class AmpCalc : public thrust::unary_function<unsigned int, devcomplex<fptype>> {
public:
    AmpCalc(unsigned int AmpIdx, unsigned int pIdx, unsigned int nPerm);
    // void setpIdx(unsigned int pIdx){_parameters = pIdx;}
    EXEC_TARGET devcomplex<fptype> operator()(thrust::tuple<int, fptype*, int> t) const;
private:
    unsigned int _nPerm;
    unsigned int _AmpIdx;
    unsigned int _parameters;
};

class NormIntegrator : public thrust::unary_function<thrust::tuple<int, int, fptype*, devcomplex<fptype>*>, fptype > {
public:
    NormIntegrator(unsigned int pIdx);
    EXEC_TARGET fptype operator()(thrust::tuple<int, int, fptype*, devcomplex<fptype>*> t) const;
private:
    unsigned int _parameters;
};


#endif

