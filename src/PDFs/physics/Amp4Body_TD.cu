/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!

TODO:
- Test lineshapes, only done for BW_DP and BW_MINT so far
- Check and implement more SF
- Currently no check if the event is even allowed in phasespace is done. This should preferably be done outside of this
class.

- Some things could be implemented differently maybe, performance should be compared for both cases.
  -For example the way Spinfactors are stored in the same array as the Lineshape values.
   Is this really worth the memory we lose by using a complex to store the SF?
*/
#include <goofit/Error.h>
#include <goofit/FitControl.h>
#include <goofit/Log.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp4Body.h>
#include <goofit/PDFs/physics/Amp4Body_TD.h>
#include <goofit/PDFs/physics/detail/Dim5.h>
#include <goofit/detail/Complex.h>
#include <mcbooster/Evaluate.h>
#include <mcbooster/EvaluateArray.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/GFunctional.h>
#include <mcbooster/GTypes.h>
#include <mcbooster/Generate.h>
#include <mcbooster/Vector4R.h>

#include <goofit/PDFs/physics/Amp4BodyGlobals.h>
#include <goofit/PDFs/physics/Amplitude.h>
#include <goofit/PDFs/physics/detail/AmpCalc_TD.h>
#include <goofit/PDFs/physics/detail/FourDblTupleAdd.h>
#include <goofit/PDFs/physics/detail/LSCalculator_TD.h>
#include <goofit/PDFs/physics/detail/NormIntegrator_TD.h>
#include <goofit/PDFs/physics/detail/NormLSCalculator_TD.h>
#include <goofit/PDFs/physics/detail/SFCalculator_TD.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

#include <cstdarg>

namespace GooFit {

struct genExp {
    fptype gamma;
    unsigned int offset;

    __host__ __device__ genExp(unsigned int c, fptype d)
        : gamma(d)
        , offset(c){};

    __host__ __device__ auto operator()(unsigned int x) const -> fptype {
        thrust::random::default_random_engine rand(1431655765);
        thrust::uniform_real_distribution<fptype> dist(0, 1);

        rand.discard(x + offset);

        return -log(dist(rand)) / gamma;
    }
};

struct genUniform {
    unsigned int offset;

    __host__ __device__ genUniform(unsigned int c)
        :offset(c){};

    __host__ __device__ auto operator()(unsigned int x) const -> fptype {
        thrust::random::default_random_engine rand(1431655765);
        thrust::uniform_real_distribution<fptype> dist(0, 3.26);

        rand.discard(x + offset);

        return dist(rand);
    }
};

struct genExp_importance_weight {
    fptype _gamma;
    __host__ __device__ genExp_importance_weight(fptype gamma)
        : _gamma(gamma){};

    __device__ auto operator()(fptype time) -> fptype {
        fptype weight = 1.0/(exp(-time * _gamma));
        //we return an importance sampling weight if our normalisation events are generated using an exponential
        return weight;
        // Should be something like: return thrust::get<1>(t) / exp(-time * gammamin);
    }
};

struct exp_functor {
    size_t tmpparam, tmpoff;
    fptype gammamin, wmax;
    exp_functor(size_t tmpparam, size_t tmpoff, fptype gammamin, fptype wmax)
        : tmpparam(tmpparam)
        , tmpoff(tmpoff)
        , gammamin(gammamin)
        , wmax(wmax) {}

    __device__ auto operator()(thrust::tuple<unsigned int, fptype, fptype *, unsigned int> t) -> fptype {
        int evtNum  = thrust::get<0>(t);
        fptype *evt = thrust::get<2>(t) + (evtNum * thrust::get<3>(t));
        // unsigned int *indices = paramIndices + tmpparam;

        // while
        //    fptype time = evt[indices[8 + indices[0]]];

        // we are grabbing the time out
        fptype time = evt[tmpparam];

        thrust::random::minstd_rand0 rand(1431655765);
        thrust::uniform_real_distribution<fptype> dist(0, 1);
        rand.discard(tmpoff + evtNum);

        return (dist(rand) * exp(-time * gammamin) * wmax) < thrust::get<1>(t);
        // Should be something like: return thrust::get<1>(t) / exp(-time * gammamin);
    }
};

struct get_pdf_val {
    __host__ __device__ fptype operator()(thrust::tuple<fptype,fptype,fptype,fptype> t){
      fptype pdf_val = thrust::get<0>(t);
      return pdf_val;
    }

  };

struct get_norm_pdf_weight{
    fptype wmax;
    __host__ __device__ get_norm_pdf_weight(fptype _wmax)
      : wmax(_wmax){};
    __host__ __device__ fptype operator()(fptype pdf_val){
      return (fptype)(pdf_val/wmax);
    }
  };

// The function of this array is to hold all the cached waves; specific
// waves are recalculated when the corresponding resonance mass or width
// changes. Note that in a multithread environment each thread needs its
// own cache, hence the '10'. Ten threads should be enough for anyone!
__device__ fpcomplex *Amps_TD[10];
/*
Constant memory array to hold specific info for amplitude calculation.
First entries are the starting points in array, necessary, because number of Lineshapes(LS) or Spinfactors(SF) can vary
|start of each Amplitude| #Linshapes | #Spinfactors | LS-indices | SF-indices|
| 1 entry per Amplitude | 1 per Amp  | 1 per Amp    | #LS in Amp| #SF in Amp|
*/
// __constant__ unsigned int AmpIndices_TD[100];

// This function gets called by the GooFit framework to get the value of the PDF.
__device__ auto device_Amp4Body_TD(fptype *evt, ParameterContainer &pc) -> fptype {
    // printf("DalitzPlot evt %i zero: %i %i %f (%f, %f).\n", evtNum, numResonances, effFunctionIdx, eff, totalAmp.real,
    // totalAmp.imag);

    int id_evt = pc.getObservable(5);

    auto evtNum = static_cast<int>(floor(0.5 + RO_CACHE(evt[id_evt])));
    // GOOFIT_TRACE("Amp4Body_TD: Number of events: {}", evtNum);

    unsigned int cacheToUse = pc.getConstant(5);
    unsigned int numAmps    = pc.getConstant(8);
    unsigned int totalSF_LS = pc.getConstant(10);

    fpcomplex AmpA(0, 0);
    fpcomplex AmpB(0, 0);
    fpcomplex amp_A, amp_B;

    int k = 0;

    for(int i = 0; i < numAmps; ++i) {
        unsigned int start = AmpIndices[i];
        unsigned int flag  = AmpIndices[start + 3 + numAmps];
        fpcomplex temp;

        // printf("flag:%i\n",flag);
        switch(flag) {
        case 0:
            amp_A = fpcomplex(pc.getParameter(4 + 2 * (i + k)), pc.getParameter(5 + 2 * (i + k)));
            temp  = Amps_TD[cacheToUse][evtNum * numAmps + i];
            AmpA += temp * amp_A;
            break;

        case 1:
            amp_B = fpcomplex(pc.getParameter(4 + 2 * (i + k)), pc.getParameter(5 + 2 * (i + k)));
            temp  = Amps_TD[cacheToUse][evtNum * numAmps + i];
            AmpB += temp * amp_B;
            break;

        case 2:
            amp_A = fpcomplex(pc.getParameter(4 + 2 * (i + k)), pc.getParameter(5 + 2 * (i + k)));
            temp  = Amps_TD[cacheToUse][evtNum * numAmps + i];
            AmpA += temp * amp_A;

            ++k;
            amp_B = fpcomplex(pc.getParameter(4 + 2 * (i + k)), pc.getParameter(5 + 2 * (i + k)));
            temp  = Amps_TD[cacheToUse][evtNum * numAmps + i];
            AmpB += temp * amp_B;
            break;
        }
    }

    int id_time  = pc.getObservable(6);
    int id_sigma = pc.getObservable(7);
    int id_wsig = pc.getObservable(8);

    fptype _tau          = pc.getParameter(0);
    fptype _xmixing      = pc.getParameter(1);
    fptype _ymixing      = pc.getParameter(2);
    fptype _SqWStoRSrate = pc.getParameter(3);
    fptype _time         = RO_CACHE(evt[id_time]);
    fptype _sigma        = RO_CACHE(evt[id_sigma]);

    //std::cout << "reading time: " << _time << "with resolution: " << _sigma << std::endl;
    //wSig is the weight associated with a per event relative efficiency
    
    fptype wSig = RO_CACHE(evt[id_wsig]);

    AmpA *= _SqWStoRSrate;
    /*printf("%i read time: %.5g x: %.5g y: %.5g \n",evtNum, _time, _xmixing, _ymixing);*/

    fptype term1    = thrust::norm(AmpA) + thrust::norm(AmpB);
    fptype term2    = thrust::norm(AmpA) - thrust::norm(AmpB);
    fpcomplex term3 = AmpA * thrust::conj(AmpB);
    // printf("%i dev %.7g %.7g %.7g %.7g\n", evtNum, norm2(AmpA), norm2(AmpB), term3.real, term3.imag);

    // increment Index
    pc.incrementIndex();

    // increment over all our lineshapes and spinfactors to get to the resolution function
    for(int i = 0; i < totalSF_LS; i++)
        pc.incrementIndex();

    // int effFunctionIdx = 12 + 2 * indices[3] + 2 * indices[4] + 2 * indices[6];
    // int resfctidx      = indices[11];
    // int resfctpar      = effFunctionIdx + 2;

    // resolution function?
    fptype ret = (*(reinterpret_cast<device_resfunction_ptr>(d_function_table[pc.funcIdx])))(
        term1, term2, term3.real(), term3.imag(), _tau, _time, _xmixing, _ymixing, _sigma, pc);

    // increment resolution function
    pc.incrementIndex();

    // efficiency function?
    fptype eff = callFunction(evt, pc);
    /*printf("%i result %.7g, eff %.7g\n",evtNum, ret, eff);*/
    ret *= wSig;

    ret *= eff;
    /*printf("in prob: %f\n", ret);*/
    return ret;
}

__device__ device_function_ptr ptr_to_Amp4Body_TD = device_Amp4Body_TD;

__host__ Amp4Body_TD::Amp4Body_TD(std::string n,
                                  std::vector<Observable> observables,
                                  DecayInfo4t decay,
                                  MixingTimeResolution *Tres,
                                  GooPdf *efficiency,
                                  Observable *mistag,
                                  unsigned int MCeventsNorm)
    : Amp4BodyBase("Amp4Body_TD", n)
    , decayInfo(decay)
    , resolution(Tres)
    , totalEventSize(observables.size() + 2) // number of observables plus eventnumber
{
    // should include m12, m34, cos12, cos34, phi, eventnumber, dtime, sigmat. In this order!
    for(auto &observable : observables) {
        registerObservable(observable);
    }

    // constantsList.push_back(decayInfo.meson_radius);

    for(double &particle_masse : decayInfo.particle_masses) {
        registerConstant(particle_masse);
    }

    static int cacheCount = 0;
    cacheToUse            = cacheCount++;
    registerParameter(decayInfo._tau);
    registerParameter(decayInfo._xmixing);
    registerParameter(decayInfo._ymixing);
    registerParameter(decayInfo._SqWStoRSrate);

    if(resolution->getDeviceFunction() < 0)
        throw GooFit::GeneralError("The resolution device function index {} must be more than 0",
                                   resolution->getDeviceFunction());

    registerConstant(cacheToUse);

    // This is the start of reading in the amplitudes and adding the lineshapes and Spinfactors to this PDF
    // This is done in this way so we don't have multiple copies of one lineshape in one pdf.
    unsigned int coeff_counter = 0;

    std::vector<Amplitude *> AmpBuffer;
    std::vector<Amplitude *> AmpsA = decayInfo.amplitudes;
    std::vector<Amplitude *> AmpsB = decayInfo.amplitudes_B;

    std::vector<unsigned int> nPermVec;
    std::vector<unsigned int> amp_idx;
    std::vector<unsigned int> amp_idx_start;

    unsigned int total_lineshapes_spinfactors = 0;

    for(auto &i : AmpsA) {
        AmpMap[i->_uniqueDecayStr] = std::make_pair(std::vector<unsigned int>(0), std::vector<unsigned int>(0));
        AmpBuffer.push_back(i);

        // adding amplitude to components
        components.push_back(i);

        // register the parameters from the amplitudes here
        registerParameter(i->_ar);
        registerParameter(i->_ai);

        auto LSvec = i->_LS;

        for(auto &LSIT : LSvec) {
            total_lineshapes_spinfactors++;
            auto found = std::find_if(
                LineShapes.begin(), LineShapes.end(), [&LSIT](const Lineshape *L) { return (*LSIT) == (*L); });

            if(found != LineShapes.end()) {
                AmpMap[i->_uniqueDecayStr].first.push_back(std::distance(LineShapes.begin(), found));
                // printf("LS %s found at %i\n",LSIT->getName().c_str(),std::distance(LineShapes.begin(), found));
            } else {
                LineShapes.push_back(LSIT);
                AmpMap[i->_uniqueDecayStr].first.push_back(LineShapes.size() - 1);
                // printf("Adding LS %s\n",LSIT->getName().c_str());
            }
        }

        auto SFvec = i->_SF;

        for(auto &SFIT : SFvec) {
            total_lineshapes_spinfactors++;
            auto found = std::find_if(
                SpinFactors.begin(), SpinFactors.end(), [&SFIT](const SpinFactor *S) { return (*SFIT) == (*S); });

            if(found != SpinFactors.end()) {
                AmpMap[i->_uniqueDecayStr].second.push_back(std::distance(SpinFactors.begin(), found));
            } else {
                SpinFactors.push_back(SFIT);
                AmpMap[i->_uniqueDecayStr].second.push_back(SpinFactors.size() - 1);
            }
        }

        unsigned int flag = 0;
        auto inB = std::find_if(AmpsB.begin(), AmpsB.end(), [AmpsA, &i](const Amplitude *A) { return *i == (*A); });

        if(inB != AmpsB.end()) {
            flag = 2;
            // printf("Found in AmpsB as well: %s\n", (*inB)->_uniqueDecayStr.c_str());
            registerParameter((*inB)->_ar);
            registerParameter((*inB)->_ai);
            ++coeff_counter;
        }

        nPermVec.push_back(i->_nPerm);
        std::vector<unsigned int> ls = AmpMap[i->_uniqueDecayStr].first;
        std::vector<unsigned int> sf = AmpMap[i->_uniqueDecayStr].second;
        amp_idx_start.push_back(amp_idx.size());
        amp_idx.push_back(ls.size());
        amp_idx.push_back(sf.size());
        amp_idx.push_back(i->_nPerm);
        amp_idx.push_back(flag);
        amp_idx.insert(amp_idx.end(), ls.begin(), ls.end());
        amp_idx.insert(amp_idx.end(), sf.begin(), sf.end());
        ++coeff_counter;
        // AmpBuffer.push_back(i);
    }

    for(auto &i : AmpsB) {
        unsigned int flag = 1;
        auto inB
            = std::find_if(AmpBuffer.begin(), AmpBuffer.end(), [AmpsB, &i](const Amplitude *A) { return *i == (*A); });

        if(inB != AmpBuffer.end())
            continue;

        AmpMap[i->_uniqueDecayStr] = std::make_pair(std::vector<unsigned int>(0), std::vector<unsigned int>(0));

        components.push_back(i);
        // AmpBuffer.push_back(i);

        registerParameter(i->_ar);
        registerParameter(i->_ai);

        auto LSvec = i->_LS;

        for(auto &LSIT : LSvec) {
            total_lineshapes_spinfactors++;
            auto found = std::find_if(
                LineShapes.begin(), LineShapes.end(), [&LSIT](const Lineshape *L) { return (*LSIT) == (*L); });

            if(found != LineShapes.end()) {
                AmpMap[i->_uniqueDecayStr].first.push_back(std::distance(LineShapes.begin(), found));
                // printf("LS %s found at %i\n",LSIT->getName().c_str(),std::distance(LineShapes.begin(), found));
            } else {
                LineShapes.push_back(LSIT);
                AmpMap[i->_uniqueDecayStr].first.push_back(LineShapes.size() - 1);
                // printf("Adding LS %s\n",LSIT->getName().c_str());
            }
        }

        auto SFvec = i->_SF;

        for(auto &SFIT : SFvec) {
            total_lineshapes_spinfactors++;
            auto found = std::find_if(
                SpinFactors.begin(), SpinFactors.end(), [&SFIT](const SpinFactor *S) { return (*SFIT) == (*S); });

            if(found != SpinFactors.end()) {
                AmpMap[i->_uniqueDecayStr].second.push_back(std::distance(SpinFactors.begin(), found));
            } else {
                SpinFactors.push_back(SFIT);
                AmpMap[i->_uniqueDecayStr].second.push_back(SpinFactors.size() - 1);
            }
        }

        nPermVec.push_back(i->_nPerm);
        ++coeff_counter;
        amp_idx_start.push_back(amp_idx.size());
        std::vector<unsigned int> ls = AmpMap[i->_uniqueDecayStr].first;
        std::vector<unsigned int> sf = AmpMap[i->_uniqueDecayStr].second;
        amp_idx.push_back(ls.size());
        amp_idx.push_back(sf.size());
        amp_idx.push_back(i->_nPerm);
        amp_idx.push_back(flag);
        amp_idx.insert(amp_idx.end(), ls.begin(), ls.end());
        amp_idx.insert(amp_idx.end(), sf.begin(), sf.end());
    }

    registerConstant(LineShapes.size());  //#LS
    registerConstant(SpinFactors.size()); //#SF
    registerConstant(components.size());  //# AMP
    registerConstant(coeff_counter); // Number of coefficients, because its not necessary to be equal to number of Amps.
    registerConstant(total_lineshapes_spinfactors);

    components.push_back(resolution);
    components.push_back(efficiency);

    if(mistag) {
        registerObservable(*mistag);
        totalEventSize = 9;
        // TODO: This needs to be registered later!
        // registerConstant(1); // Flags existence of mistag
    }

    // In case the resolution function needs parameters, this registers them.
    // resolution->createParameters(pindices, this);

    registerFunction("ptr_to_Amp4Body_TD", ptr_to_Amp4Body_TD);

    initialize();

    //Integrator   = new NormIntegrator_TD();
    fprintf(stderr,"Special integral: %d",Amp4Body_TD::specialIntegral);
    Integrator = new NormIntegrator_TD(Amp4Body_TD::specialIntegral);
    redoIntegral = new bool[LineShapes.size()];
    cachedMasses = new fptype[LineShapes.size()];
    cachedWidths = new fptype[LineShapes.size()];

    for(int i = 0; i < LineShapes.size(); i++) {
        lscalculators.push_back(new LSCalculator_TD());
    }

    MEMCPY_TO_SYMBOL(
        AmpIndices, &(amp_idx_start[0]), amp_idx_start.size() * sizeof(unsigned int), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(AmpIndices,
                     &(amp_idx[0]),
                     amp_idx.size() * sizeof(unsigned int),
                     amp_idx_start.size() * sizeof(unsigned int),
                     cudaMemcpyHostToDevice);

    for(int i = 0; i < SpinFactors.size(); ++i) {
        sfcalculators.push_back(new SFCalculator_TD());
    }

    for(int i = 0; i < components.size() - 2; ++i) {
        AmpCalcs.push_back(new AmpCalc_TD(nPermVec[i], amp_idx_start[i]));
    }

    // fprintf(stderr,"#Amp's %i, #LS %i, #SF %i \n", AmpMap.size(), components.size()-1, SpinFactors.size() );

    std::vector<mcbooster::GReal_t> masses(decayInfo.particle_masses.begin() + 1, decayInfo.particle_masses.end());
    mcbooster::PhaseSpace phsp(decayInfo.particle_masses[0], masses, MCeventsNorm, generation_offset);
    phsp.Generate(mcbooster::Vector4R(decayInfo.particle_masses[0], 0.0, 0.0, 0.0));
    phsp.Unweight();

    auto nAcc                     = phsp.GetNAccepted();
    mcbooster::BoolVector_d flags = phsp.GetAccRejFlags();
    auto d1                       = phsp.GetDaughters(0);
    auto d2                       = phsp.GetDaughters(1);
    auto d3                       = phsp.GetDaughters(2);
    auto d4                       = phsp.GetDaughters(3);

    auto zip_begin = thrust::make_zip_iterator(thrust::make_tuple(d1.begin(), d2.begin(), d3.begin(), d4.begin()));
    auto zip_end   = zip_begin + d1.size();
    auto new_end   = thrust::remove_if(zip_begin, zip_end, flags.begin(), thrust::logical_not<bool>());

    d1.erase(thrust::get<0>(new_end.get_iterator_tuple()), d1.end());
    d2.erase(thrust::get<1>(new_end.get_iterator_tuple()), d2.end());
    d3.erase(thrust::get<2>(new_end.get_iterator_tuple()), d3.end());
    d4.erase(thrust::get<3>(new_end.get_iterator_tuple()), d4.end());

    mcbooster::ParticlesSet_d pset(4);
    pset[0] = &d1;
    pset[1] = &d2;
    pset[2] = &d3;
    pset[3] = &d4;

    norm_M12        = mcbooster::RealVector_d(nAcc);
    norm_M34        = mcbooster::RealVector_d(nAcc);
    norm_CosTheta12 = mcbooster::RealVector_d(nAcc);
    norm_CosTheta34 = mcbooster::RealVector_d(nAcc);
    norm_phi        = mcbooster::RealVector_d(nAcc);

    //normalistion generated decay time
    norm_dtime = mcbooster::RealVector_d(nAcc);
    //efficiency weight from BDT for the normalisation events
    norm_eff = mcbooster::RealVector_d(nAcc);
    //weight associated from the pdf value in place of an accept-reject
    //norm_pdf_weight = mcbooster::RealVector_d(nAcc);
    //weight from importance sampling the decay time distribution
    norm_importance_weight = mcbooster::RealVector_d(nAcc);
    //Do this straight after intialisation
    thrust::counting_iterator<unsigned int> index_sequence_begin(0);

    fptype tau      = parametersList[0].getValue();
    //adding the x mixing term to see if this affects the accept reject numbers in gammamin
    //fptype xmixing = parametersList[1].getValue();
    fptype ymixing  = parametersList[2].getValue();
    fptype gammamin = 1.0 / tau - fabs(ymixing) / tau;
    //fill the normalisation decay times with zero initially to calculate the value of the PDF at t=0
    //thrust::fill(norm_dtime.begin(),norm_dtime.end(),0.);
    //fill the normalisation weights with 1
    //thrust::fill(norm_pdf_weight.begin(),norm_pdf_weight.end(),1.);
    thrust::fill(norm_importance_weight.begin(),norm_importance_weight.end(),1.);
    //fill the normalisation vectors with ones to avoid any errors involving multiplying by 0
    thrust::fill(norm_eff.begin(),norm_eff.end(),1.);

    //generate uniformly distrubted decay times if not using importance sampling
    //thrust::transform(index_sequence_begin, index_sequence_begin + nAcc, norm_dtime.begin(), genUniform(generation_offset));
    
    thrust::transform(index_sequence_begin, index_sequence_begin + nAcc, norm_dtime.begin(), genExp(generation_offset, gammamin));
    thrust::transform(norm_dtime.begin(), norm_dtime.end(), norm_importance_weight.begin(), genExp_importance_weight(gammamin));
    
    mcbooster::RealVector_h norm_importance_weight_h(norm_importance_weight);
    for(int i = 0; i < 10;i++){
        fprintf(stderr,"importance weight %.3g\n",norm_importance_weight_h[i]);
    }
    
    mcbooster::VariableSet_d VarSet(5);
    VarSet[0] = &norm_M12;
    VarSet[1] = &norm_M34;
    VarSet[2] = &norm_CosTheta12;
    VarSet[3] = &norm_CosTheta34;
    VarSet[4] = &norm_phi;

    Dim5 eval = Dim5();
    mcbooster::EvaluateArray<Dim5>(eval, pset, VarSet);

    norm_SF  = mcbooster::RealVector_d(nAcc * SpinFactors.size());
    norm_LS  = mcbooster::mc_device_vector<fpcomplex>(nAcc * (components.size() - 1));
    MCevents = nAcc;

    setSeparateNorm();
}

__host__ void Amp4Body_TD::populateArrays() {
    PdfBase::populateArrays();

    // TODO: We need to expand populateArrays so we handle components correctly!
    efficiencyFunction = host_function_table.size() - 1;
}

// makes the arrays to chache the lineshape values and spinfactors in CachedResSF and the values of the amplitudes in
// cachedAMPs
// I made the choice to have spinfactors necxt to the values of the lineshape in memory. I waste memory by doing this
// because a spinfactor is saved as complex
// It would be nice to test if this is better than having the spinfactors stored seperately.
__host__ void Amp4Body_TD::setDataSize(unsigned int dataSize, unsigned int evtSize) {
    // Default 3 is m12, m13, evtNum for DP 2dim, 4-body decay has 5 independent vars plus evtNum = 6
    totalEventSize = evtSize;
    if(totalEventSize < 3)
        throw GooFit::GeneralError("totalEventSize {} must be 3 or more", totalEventSize);

    if(cachedResSF)
        delete cachedResSF;

    if(cachedAMPs)
        delete cachedAMPs;

    numEntries  = dataSize;
    cachedResSF = new thrust::device_vector<fpcomplex>(dataSize * (LineShapes.size() + SpinFactors.size()));
    void *dummy = thrust::raw_pointer_cast(cachedResSF->data());
    MEMCPY_TO_SYMBOL(cResSF_TD, &dummy, sizeof(fpcomplex *), cacheToUse * sizeof(fpcomplex *), cudaMemcpyHostToDevice);

    GOOFIT_TRACE("cachedAmps size:{}, {}, {}, {}, {}",
                 dataSize,
                 AmpCalcs.size(),
                 components.size() - 1,
                 SpinFactors.size(),
                 LineShapes.size());

    cachedAMPs   = new thrust::device_vector<fpcomplex>(dataSize * (AmpCalcs.size()));
    void *dummy2 = thrust::raw_pointer_cast(cachedAMPs->data());
    MEMCPY_TO_SYMBOL(Amps_TD, &dummy2, sizeof(fpcomplex *), cacheToUse * sizeof(fpcomplex *), cudaMemcpyHostToDevice);

    setForceIntegrals();
}

// this is where the actual magic happens. This function does all the calculations!
__host__ auto Amp4Body_TD::normalize() -> fptype {
    if(cachedResSF == nullptr)
        throw GeneralError("You must call dp.setDataSize(currData.getNumEvents(), N) first!");
    // fprintf(stderr, "start normalize\n");
    recursiveSetNormalization(1.0); // Not going to normalize efficiency,
    // so set normalization factor to 1 so it doesn't get multiplied by zero.
    // Copy at this time to ensure that the SpecialResonanceCalculators, which need the efficiency,
    // don't get zeroes through multiplying by the normFactor.
    host_normalizations.sync(d_normalizations);

    // check if MINUIT changed any parameters and if so remember that so we know
    // we need to recalculate that lineshape and every amp, that uses that lineshape
    for(unsigned int i = 0; i < components.size() - 1; ++i) {
        redoIntegral[i] = forceRedoIntegrals;

        if(!(components[i]->parametersChanged()))
            continue;

        redoIntegral[i] = true;
    }

    SpinsCalculated    = !forceRedoIntegrals;
    forceRedoIntegrals = false;

    // just some thrust iterators for the calculation.
    thrust::constant_iterator<fptype *> dataArray(dev_event_array);
    thrust::constant_iterator<int> eventSize(totalEventSize);
    thrust::counting_iterator<int> eventIndex(0);

    // Calculate spinfactors only once for normalization events and real events
    // strided_range is a template implemented in DalitsPlotHelpers.hh
    // it basically goes through the array by increasing the pointer by a certain amount instead of just one step.
    if(!SpinsCalculated) {
        for(int i = 0; i < SpinFactors.size(); ++i) {
            unsigned int offset = LineShapes.size();
            unsigned int stride = LineShapes.size() + SpinFactors.size();

            GOOFIT_TRACE("SpinFactors - stride: {}", stride);
            sfcalculators[i]->setDalitzId(getFunctionIndex());
            sfcalculators[i]->setSpinFactorId(SpinFactors[i]->getFunctionIndex());

            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, dataArray, eventSize)),
                strided_range<thrust::device_vector<fpcomplex>::iterator>(
                    cachedResSF->begin() + offset + i, cachedResSF->end(), stride)
                    .begin(),
                *(sfcalculators[i]));

            if(!generation_no_norm) {
                NormSpinCalculator_TD nsc = NormSpinCalculator_TD();

                nsc.setDalitzId(getFunctionIndex());
                nsc.setSpinFactorId(SpinFactors[i]->getFunctionIndex());

                thrust::transform(
                    thrust::make_zip_iterator(thrust::make_tuple(norm_M12.begin(),
                                                                 norm_M34.begin(),
                                                                 norm_CosTheta12.begin(),
                                                                 norm_CosTheta34.begin(),
                                                                 norm_phi.begin())),
                    thrust::make_zip_iterator(thrust::make_tuple(
                        norm_M12.end(), norm_M34.end(), norm_CosTheta12.end(), norm_CosTheta34.end(), norm_phi.end())),
                    (norm_SF.begin() + (i * MCevents)),
                    nsc);
            }
        }

        SpinsCalculated = true;
    }

    // fprintf(stderr, "normalize after spins\n");

    // this calculates the values of the lineshapes and stores them in the array. It is recalculated every time
    // parameters change.
    for(int i = 0; i < LineShapes.size(); ++i) {
        if(redoIntegral[i]) {
            lscalculators[i]->setDalitzId(getFunctionIndex());
            lscalculators[i]->setResonanceId(LineShapes[i]->getFunctionIndex());

            unsigned int stride = LineShapes.size() + SpinFactors.size();

            GOOFIT_TRACE("LineShape[{}] - stride: {}", LineShapes[i]->getName().c_str(), stride);
            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, dataArray, eventSize)),
                strided_range<thrust::device_vector<fpcomplex>::iterator>(
                    cachedResSF->begin() + i, cachedResSF->end(), stride)
                    .begin(),
                *(lscalculators[i]));
        }
    }

    // fprintf(stderr, "normalize after LS\n");

    // this is a little messy but it basically checks if the amplitude includes one of the recalculated lineshapes and
    // if so recalculates that amplitude
    // auto AmpMapIt = AmpMap.begin();

    for(int i = 0; i < components.size() - 2; ++i) {
        bool redo = false;
        for(unsigned int j = 0; j < components.size() - 2; j++) {
            if(!redoIntegral[j])
                continue;
            redo = true;
            break;
        }

        if(redo) {
            AmpCalcs[i]->setDalitzId(getFunctionIndex());

            thrust::transform(eventIndex,
                              eventIndex + numEntries,
                              strided_range<thrust::device_vector<fpcomplex>::iterator>(
                                  cachedAMPs->begin() + i, cachedAMPs->end(), AmpCalcs.size())
                                  .begin(),
                              *(AmpCalcs[i]));
        }
    }

    // fprintf(stderr, "normalize after Amps\n");

    // lineshape value calculation for the normalization, also recalculated every time parameter change
    if(!generation_no_norm) {
        for(int i = 0; i < LineShapes.size(); ++i) {
            if(!redoIntegral[i])
                continue;

            NormLSCalculator_TD ns;
            ns.setDalitzId(getFunctionIndex());
            ns.setResonanceId(LineShapes[i]->getFunctionIndex());

            thrust::transform(
                thrust::make_zip_iterator(thrust::make_tuple(norm_M12.begin(),
                                                             norm_M34.begin(),
                                                             norm_CosTheta12.begin(),
                                                             norm_CosTheta34.begin(),
                                                             norm_phi.begin())),
                thrust::make_zip_iterator(thrust::make_tuple(
                    norm_M12.end(), norm_M34.end(), norm_CosTheta12.end(), norm_CosTheta34.end(), norm_phi.end())),
                (norm_LS.begin() + (i * MCevents)),
                ns);
        }
    }

    thrust::constant_iterator<fptype *> normSFaddress(thrust::raw_pointer_cast(norm_SF.data()));
    thrust::constant_iterator<fpcomplex *> normLSaddress(thrust::raw_pointer_cast(norm_LS.data()));
    thrust::constant_iterator<int> NumNormEvents(MCevents);

    // this does the rest of the integration with the cached lineshape and spinfactor values for the normalization
    // events
    auto ret = 1.0;

    if(!generation_no_norm) {
        thrust::tuple<fptype, fptype, fptype, fptype> dummy(0, 0, 0, 0);
        FourDblTupleAdd MyFourDoubleTupleAdditionFunctor;
        thrust::tuple<fptype, fptype, fptype, fptype> sumIntegral;

        Integrator->setDalitzId(getFunctionIndex());

        fprintf(stderr,"Evaluating numerical normalization\n");
        sumIntegral = thrust::transform_reduce(
            thrust::make_zip_iterator(thrust::make_tuple(eventIndex, NumNormEvents, normSFaddress, normLSaddress,
        norm_dtime.begin(),norm_eff.begin(),norm_importance_weight.begin())),
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex + MCevents, NumNormEvents, normSFaddress, normLSaddress,
        norm_dtime.end(),norm_eff.end(),norm_importance_weight.end())),
        *Integrator, dummy, MyFourDoubleTupleAdditionFunctor);  
        fprintf(stderr,"Finished evaluating numerical normalization\n");
    // GOOFIT_TRACE("sumIntegral={}", sumIntegral);

        // printf("normalize A2/#evts , B2/#evts: %.5g, %.5g\n",thrust::get<0>(sumIntegral)/MCevents,
        // thrust::get<1>(sumIntegral)/MCevents);
        fptype tau     = parametersList[0];
        fptype xmixing = parametersList[1];
        fptype ymixing = parametersList[2];
        if(specialIntegral){
            const double uniformNorm = 1.;
            //const double uniformNorm = 3.26 - 0.18;
            ret = thrust::get<0>(sumIntegral) * uniformNorm;
            auto ratioNormBA = thrust::get<2>(sumIntegral) / thrust::get<1>(sumIntegral);
            //fprintf(stderr, "SpecInt normalize A2/#evts , B2/#evts, ratio: %.5g, %.5g, %.5g\n",
            //      thrust::get<1>(sumIntegral)/MCevents, thrust::get<2>(sumIntegral)/MCevents, ratioNormBA);
          }
          else {
            ret = resolution->normalization(thrust::get<0>(sumIntegral),
            thrust::get<1>(sumIntegral),
            thrust::get<2>(sumIntegral),
            thrust::get<3>(sumIntegral),
            tau,
            xmixing,
            ymixing);

          }
 

        // MCevents is the number of normalization events.
        ret /= MCevents;
        if(specialIntegral){
            printf("normalization value:%.7g \n",ret);
          }
    }

    host_normalizations[normalIdx + 1] = 1.0 / ret;
    cachedNormalization                = 1.0 / ret;
    // printf("end of normalize %f\n", ret);
    return ret;
}

__host__ auto Amp4Body_TD::GenerateSig(unsigned int numEvents, int seed) -> std::
    tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::BoolVector_h> {
    initialize();
    copyParams();

    std::vector<mcbooster::GReal_t> masses(decayInfo.particle_masses.begin() + 1, decayInfo.particle_masses.end());
    mcbooster::PhaseSpace phsp(decayInfo.particle_masses[0], masses, numEvents, generation_offset);
    if(seed != 0)
        phsp.SetSeed(seed);
    else
        GOOFIT_INFO("Current generator seed {}, offset {}", phsp.GetSeed(), generation_offset);

    phsp.Generate(mcbooster::Vector4R(decayInfo.particle_masses[0], 0.0, 0.0, 0.0));

    phsp.Unweight();

    auto nAcc                     = phsp.GetNAccepted();
    mcbooster::BoolVector_d flags = phsp.GetAccRejFlags();
    auto d1                       = phsp.GetDaughters(0);
    auto d2                       = phsp.GetDaughters(1);
    auto d3                       = phsp.GetDaughters(2);
    auto d4                       = phsp.GetDaughters(3);

    auto zip_begin = thrust::make_zip_iterator(thrust::make_tuple(d1.begin(), d2.begin(), d3.begin(), d4.begin()));
    auto zip_end   = zip_begin + d1.size();
    auto new_end   = thrust::remove_if(zip_begin, zip_end, flags.begin(), thrust::logical_not<bool>());

    flags = mcbooster::BoolVector_d();

    d1.erase(thrust::get<0>(new_end.get_iterator_tuple()), d1.end());
    d2.erase(thrust::get<1>(new_end.get_iterator_tuple()), d2.end());
    d3.erase(thrust::get<2>(new_end.get_iterator_tuple()), d3.end());
    d4.erase(thrust::get<3>(new_end.get_iterator_tuple()), d4.end());

    mcbooster::ParticlesSet_d pset(4);
    pset[0] = &d1;
    pset[1] = &d2;
    pset[2] = &d3;
    pset[3] = &d4;

    auto SigGen_M12_d        = mcbooster::RealVector_d(nAcc);
    auto SigGen_M34_d        = mcbooster::RealVector_d(nAcc);
    auto SigGen_CosTheta12_d = mcbooster::RealVector_d(nAcc);
    auto SigGen_CosTheta34_d = mcbooster::RealVector_d(nAcc);
    auto SigGen_phi_d        = mcbooster::RealVector_d(nAcc);
    auto dtime_d             = mcbooster::RealVector_d(nAcc);

    thrust::counting_iterator<unsigned int> index_sequence_begin(0);

    fptype tau      = parametersList[0].getValue();
    fptype ymixing  = parametersList[2].getValue();
    fptype gammamin = 1.0 / tau - fabs(ymixing) / tau;

    thrust::transform(
        index_sequence_begin, index_sequence_begin + nAcc, dtime_d.begin(), genExp(generation_offset, gammamin));

    mcbooster::VariableSet_d VarSet_d(5);
    VarSet_d[0] = &SigGen_M12_d;
    VarSet_d[1] = &SigGen_M34_d;
    VarSet_d[2] = &SigGen_CosTheta12_d;
    VarSet_d[3] = &SigGen_CosTheta34_d;
    VarSet_d[4] = &SigGen_phi_d;

    Dim5 eval = Dim5();
    mcbooster::EvaluateArray<Dim5>(eval, pset, VarSet_d);

    auto h1 = new mcbooster::Particles_h(d1);
    auto h2 = new mcbooster::Particles_h(d2);
    auto h3 = new mcbooster::Particles_h(d3);
    auto h4 = new mcbooster::Particles_h(d4);

    mcbooster::ParticlesSet_h ParSet(4);
    ParSet[0] = h1;
    ParSet[1] = h2;
    ParSet[2] = h3;
    ParSet[3] = h4;

    auto SigGen_M12_h        = new mcbooster::RealVector_h(SigGen_M12_d);
    auto SigGen_M34_h        = new mcbooster::RealVector_h(SigGen_M34_d);
    auto SigGen_CosTheta12_h = new mcbooster::RealVector_h(SigGen_CosTheta12_d);
    auto SigGen_CosTheta34_h = new mcbooster::RealVector_h(SigGen_CosTheta34_d);
    auto SigGen_phi_h        = new mcbooster::RealVector_h(SigGen_phi_d);
    auto dtime_h             = new mcbooster::RealVector_h(dtime_d);

    mcbooster::VariableSet_h VarSet(6);
    VarSet[0] = SigGen_M12_h;
    VarSet[1] = SigGen_M34_h;
    VarSet[2] = SigGen_CosTheta12_h;
    VarSet[3] = SigGen_CosTheta34_h;
    VarSet[4] = SigGen_phi_h;
    VarSet[5] = dtime_h;

    phsp.FreeResources();

    auto DS = new mcbooster::RealVector_d(8 * nAcc);
    thrust::counting_iterator<int> eventNumber(0);

#pragma unroll

    for(int i = 0; i < 5; ++i) {
        mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr(DS->begin() + i, DS->end(), 8);
        thrust::copy(VarSet_d[i]->begin(), VarSet_d[i]->end(), sr.begin());
    }

    mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr(DS->begin() + 5, DS->end(), 8);
    thrust::copy(eventNumber, eventNumber + nAcc, sr.begin());

    mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr2(DS->begin() + 6, DS->end(), 8);
    thrust::fill_n(sr2.begin(), nAcc, 0);

    dev_event_array = thrust::raw_pointer_cast(DS->data());
    setDataSize(nAcc, 8);

    generation_no_norm = true; // we need no normalization for generation, but we do need to make sure that norm = 1;
    SigGenSetIndices();
    normalize();
    setForceIntegrals();
    host_normalizations.sync(d_normalizations);

    thrust::device_vector<fptype> weights(nAcc);
    thrust::constant_iterator<int> eventSize(8);
    thrust::constant_iterator<fptype *> arrayAddress(dev_event_array);
    thrust::counting_iterator<int> eventIndex(0);

    auto fc = fitControl;
    setFitControl(std::make_shared<ProbFit>());
    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
                      thrust::make_zip_iterator(thrust::make_tuple(eventIndex + nAcc, arrayAddress, eventSize)),
                      weights.begin(),
                      *logger);
    cudaDeviceSynchronize();

    fptype wmax = 1.1 * (fptype)*thrust::max_element(weights.begin(), weights.end());

    if(wmax > maxWeight && maxWeight != 0)
        fprintf(stderr,
                "WARNING: you just encountered a higher maximum weight than observed in previous iterations.\n"
                "WARNING: Consider recalculating your AccRej flags and acceping based upon these.\n"
                "WARNING: previous weight: %.4g, new weight: %.4g\n",
                maxWeight,
                wmax);

    maxWeight = wmax > maxWeight ? wmax : maxWeight;

    thrust::copy(dtime_d.begin(), dtime_d.end(), sr2.begin());

    dtime_d = mcbooster::RealVector_d();
    thrust::device_vector<fptype> results(nAcc);

    thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
                      thrust::make_zip_iterator(thrust::make_tuple(eventIndex + nAcc, arrayAddress, eventSize)),
                      results.begin(),
                      *logger);

    setFitControl(fc);

    cudaDeviceSynchronize();

    thrust::device_vector<fptype> flag2(nAcc); // Should be weight2
    thrust::counting_iterator<mcbooster::GLong_t> first(0);
    // thrust::counting_iterator<mcbooster::GLong_t> last = first + nAcc;

    // we do not want to copy the whole class to the GPU so capturing *this is not a great option
    // therefore perpare local copies to capture the variables we need
    unsigned int tmpoff   = generation_offset;
    unsigned int tmpparam = 6;
    wmax                  = maxWeight;

    thrust::transform(
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex, results.begin(), arrayAddress, eventSize)),
        thrust::make_zip_iterator(thrust::make_tuple(eventIndex + nAcc, results.end(), arrayAddress, eventSize)),
        flag2.begin(), // Should be weight2
        exp_functor(tmpparam, tmpoff, gammamin, wmax));

    // set weights with weight2

    // Maybe max goes here now

    // Include code from Amp4Body regular here

    gooFree(dev_event_array);

    auto weights_h = mcbooster::RealVector_h(weights);
    auto results_h = mcbooster::RealVector_h(results);
    auto flags_h   = mcbooster::BoolVector_h(flag2);
    cudaDeviceSynchronize();

    return std::make_tuple(ParSet, VarSet, weights_h, flags_h);
}

} // namespace GooFit
