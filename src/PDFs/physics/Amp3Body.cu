#include <mcbooster/Evaluate.h>
#include <mcbooster/EvaluateArray.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/GFunctional.h>
#include <mcbooster/GTypes.h>
#include <mcbooster/Generate.h>
#include <mcbooster/Vector4R.h>

#include <goofit/Error.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/Amp3Body.h>
#include <goofit/PDFs/physics/Amp3BodyBase.h>
#include <goofit/PDFs/physics/detail/Dim2.h>
#include <goofit/PDFs/physics/detail/SpecialResonanceCalculator.h>
#include <goofit/PDFs/physics/detail/SpecialResonanceIntegrator.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>
#include <goofit/PDFs/physics/resonances/ResonanceUtils.h>
#include <goofit/detail/Complex.h>

#include <thrust/copy.h>
#include <thrust/transform_reduce.h>

#include <array>
#include <vector>
#include <fstream>

namespace GooFit {

// Functor used for fit fraction sum
struct CoefSumFunctor {
    fpcomplex coef_i;
    fpcomplex coef_j;

    CoefSumFunctor(fpcomplex coef_i, fpcomplex coef_j)
        : coef_i(coef_i)
        , coef_j(coef_j) {}

    __device__ auto operator()(thrust::tuple<fpcomplex, fpcomplex> val) -> fptype {
        return (coef_i * thrust::conj<fptype>(coef_j) * thrust::get<0>(val) * thrust::conj<fptype>(thrust::get<1>(val)))
            .real();
    }
};

constexpr int resonanceOffset_DP = 4; // Offset of the first resonance into the parameter index array
// Offset is number of parameters, constant index, number of resonances (not calculable
// from nP because we don't know what the efficiency might need), and cache index. Efficiency
// parameters are after the resonance information.

// The function of this array is to hold all the cached waves; specific
// waves are recalculated when the corresponding resonance mass or width
// changes. Note that in a multithread environment each thread needs its
// own cache, hence the '10'. Ten threads should be enough for anyone!

// NOTE: This is does not support ten instances (ten threads) of resoncances now, only one set of resonances.
// this needs to be large enough to hold all samples
constexpr int NUMRES = 20;


__device__ fpcomplex *cResonances[NUMRES * 20];

__device__ inline auto parIndexFromResIndex_DP(int resIndex) -> int {
    return resonanceOffset_DP + resIndex * resonanceSize;
}

__device__ auto device_DalitzPlot(fptype *evt, ParameterContainer &pc) -> fptype {
    int num_obs = pc.getNumObservables();
    int id_m13  = pc.getObservable(0);
    int id_m23  = pc.getObservable(1);
    int id_num  = pc.getObservable(2);

    fptype m13 = evt[id_m13];
    fptype m23 = evt[id_m23];

    unsigned int numResonances = pc.getConstant(0);
    unsigned int cacheToUse    = pc.getConstant(1);

    if(!inDalitz2(m13, m23, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass)) {
        pc.incrementIndex(1, numResonances * 2, 2, num_obs, 1);

        // loop over resonances and efficiency functions
        for(int i = 0; i < numResonances; i++)
            pc.incrementIndex();
        // increment the efficiency function
        pc.incrementIndex();
        return 0;
    }

    // if(c_SymDp && m13<m23)
    //     return 0;

    fptype evtIndex = evt[id_num];

    auto evtNum = static_cast<int>(floor(0.5 + evtIndex));

    fpcomplex totalAmp(0, 0);

    for(int i = 0; i < numResonances; ++i) {
        fptype mag = pc.getParameter(i * 2);
        fptype phs = pc.getParameter(i * 2 + 1);

        fpcomplex amp = thrust::polar(mag, phs);
        fpcomplex me = RO_CACHE(cResonances[i + (NUMRES * cacheToUse)][evtNum]);

        totalAmp += amp * me;
    }

    fptype ret = thrust::norm(totalAmp);
    pc.incrementIndex(1, numResonances * 2, 2, num_obs, 1);

    // loop to efficiency idx
    for(int i = 0; i < numResonances; i++)
        pc.incrementIndex();

    fptype eff = callFunction(evt, pc);
    ret *= eff;

    return ret;
}

int Amp3Body::cacheCount                         = 0;
__device__ device_function_ptr ptr_to_DalitzPlot = device_DalitzPlot;

__host__ Amp3Body::Amp3Body(
    std::string n, Observable m13, Observable m23, EventNumber eventNumber, DecayInfo3 decay, GooPdf *efficiency)
    : Amp3BodyBase("Amp3Body", n, m13, m23, eventNumber)
    , decayInfo(decay)
    , _m13(m13)
    , _m23(m23)
    , _eventNumber(eventNumber)
    , dalitzNormRange(nullptr)
    , forceRedoIntegrals(true)
    , totalEventSize(3) // Default 3 = m12, m13, evtNum
    , cacheToUse(0)
    , integrators(nullptr)
    , SymDp(&decay.SymDp)
    , calculators(nullptr) {
    for(auto &cachedWave : cachedWaves)
        cachedWave = nullptr;

    // Passing values to the defined constants.  Rather than push into list, which means each resonance
    MEMCPY_TO_SYMBOL(c_motherMass, &decay.motherMass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug1Mass, &decay.daug1Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug2Mass, &decay.daug2Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_daug3Mass, &decay.daug3Mass, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_meson_radius, &decay.meson_radius, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_mother_meson_radius, &decay.mother_meson_radius, sizeof(fptype), 0, cudaMemcpyHostToDevice);
    MEMCPY_TO_SYMBOL(c_SymDp, &decay.SymDp, sizeof(bool), 0, cudaMemcpyHostToDevice);

    // registered to 0 position
    registerConstant(decayInfo.resonances.size());

    // TODO increase after registerConstant?
    cacheToUse = cacheCount++;
    // registered to 1 position
    registerConstant(cacheToUse);

    for(auto &resonance : decayInfo.resonances) {
        

        // registering 2 parameters
        registerParameter(resonance->amp_real);
        registerParameter(resonance->amp_imag);
        components.push_back(resonance);
    }

    components.push_back(efficiency);

    registerFunction("ptr_to_DalitzPlot", ptr_to_DalitzPlot);

    initialize();

    redoIntegral   = new bool[decayInfo.resonances.size()];
    cachedMasses   = new fptype[decayInfo.resonances.size()];
    cachedWidths   = new fptype[decayInfo.resonances.size()];
    integrators    = new SpecialResonanceIntegrator **[decayInfo.resonances.size()];
    calculators    = new SpecialResonanceCalculator *[decayInfo.resonances.size()];

    for(int i = 0; i < decayInfo.resonances.size(); ++i) {
        redoIntegral[i] = true;
        cachedMasses[i] = -1;
        cachedWidths[i] = -1;
        integrators[i]  = new SpecialResonanceIntegrator *[decayInfo.resonances.size()];
        calculators[i]  = new SpecialResonanceCalculator(parameters, i);

        for(int j = 0; j < decayInfo.resonances.size(); ++j) {
            integrals.push_back(thrust::make_tuple(fpcomplex(0.,0.),fpcomplex(0.,0.)));
            integrators[i][j]    = new SpecialResonanceIntegrator(parameters, i, j);
        }
    }

    setSeparateNorm();
}

void Amp3Body::populateArrays() {
    PdfBase::populateArrays();

    // save our efficiency function.  Resonance's are saved first, then the efficiency function.  Take -1 as efficiency!
    efficiencyFunction = host_function_table.size() - 1;
}
__host__ void Amp3Body::setDataSize(unsigned int dataSize, unsigned int evtSize, unsigned int offset) {
    // Default 3 is m12, m13, evtNum
    totalEventSize = evtSize;
    if(totalEventSize < 3)
        throw GooFit::GeneralError("totalEventSize {} must be 3 or more", totalEventSize);

    // if (cachedWaves) delete cachedWaves;
    if(cachedWaves[0]) {
        for(auto &cachedWave : cachedWaves) {
            delete cachedWave;
            cachedWave = nullptr;
        }
    }

    numEntries  = dataSize;
    eventOffset = offset;

    for(int i = 0; i < NUMRES; i++) {
#ifdef GOOFIT_MPI
        cachedWaves[i] = new thrust::device_vector<fpcomplex>(m_iEventsPerTask);
#else
        cachedWaves[i] = new thrust::device_vector<fpcomplex>(dataSize);
#endif
        void *dummy = thrust::raw_pointer_cast(cachedWaves[i]->data());
        MEMCPY_TO_SYMBOL(cResonances,
                         &dummy,
                         sizeof(fpcomplex *),
                         ((NUMRES * cacheToUse) + i) * sizeof(fpcomplex *),
                         cudaMemcpyHostToDevice);
    }

    setForceIntegrals();
}

__host__ auto Amp3Body::normalize() -> fptype {
    recursiveSetNormalization(1.0); // Not going to normalize efficiency,
    // so set normalization factor to 1 so it doesn't get multiplied by zero.
    // Copy at this time to ensure that the SpecialResonanceCalculators, which need the efficiency,
    // don't get zeroes through multiplying by the normFactor.
    // we need to update the normal here, as values are used at this point.

    host_normalizations.sync(d_normalizations);

    int totalBins = _m13.getNumBins() * _m23.getNumBins();
    fptype binSizeFactor = _m13.getBinSize() * _m23.getBinSize();

    if(!dalitzNormRange) {
        gooMalloc((void **)&dalitzNormRange, 6 * sizeof(fptype));
    }

    // This line runs once
    static std::array<fptype, 6> host_norms{{0, 0, 0, 0, 0, 0}};

    std::array<fptype, 6> current_host_norms{{_m13.getLowerLimit(),
                                              _m13.getUpperLimit(),
                                              static_cast<fptype>(_m13.getNumBins()),
                                              _m23.getLowerLimit(),
                                              _m23.getUpperLimit(),
                                              static_cast<fptype>(_m23.getNumBins())}};

    if(host_norms != current_host_norms) {
        host_norms = current_host_norms;
    }
    MEMCPY(dalitzNormRange, host_norms.data(), 6 * sizeof(fptype), cudaMemcpyHostToDevice);
    for(unsigned int i = 0; i < decayInfo.resonances.size(); ++i) {
        redoIntegral[i] = forceRedoIntegrals;

        if(!(decayInfo.resonances[i]->parametersChanged()))
            continue;

        redoIntegral[i] = true;
    }

    forceRedoIntegrals = false;

    // Only do this bit if masses or widths have changed.
    thrust::constant_iterator<fptype *> arrayAddress(dalitzNormRange);
    thrust::counting_iterator<int> binIndex(0);

    // NB, SpecialResonanceCalculator assumes that fit is unbinned!
    // And it needs to know the total event size, not just observables
    // for this particular PDF component.
    thrust::constant_iterator<fptype *> dataArray(dev_event_array);
    thrust::constant_iterator<int> eventSize(totalEventSize);
    thrust::counting_iterator<int> eventIndex(eventOffset);

    const size_t n_res= decayInfo.resonances.size();

    auto reduce_func = [] __host__ __device__ (const thrust::tuple<fpcomplex, fpcomplex>& x, const thrust::tuple<fpcomplex, fpcomplex>& y)
    {
        return thrust::make_tuple(thrust::get<0>(x) + thrust::get<0>(y), thrust::get<1>(x) + thrust::get<1>(y));
    };

    for(int i = 0; i < n_res; ++i) {
        for(int j = 0; j < n_res; ++j) {
            if((!redoIntegral[i]) && (!redoIntegral[j])){
                continue;
            }
                //std::cout << "redoIntegral " << decayInfo.resonances[i]->getName() << " " << decayInfo.resonances[j]->getName()<< " " << redoIntegral[i] << " " << redoIntegral[j] << " index = " << i*n_res + j << "\n";
                integrators[i][j]->setDalitzIndex(getFunctionIndex());
                integrators[i][j]->setResonanceIndex(decayInfo.resonances[i]->getFunctionIndex());
                integrators[i][j]->setEfficiencyIndex(decayInfo.resonances[j]->getFunctionIndex());
                thrust::constant_iterator<int> effFunc(efficiencyFunction);
                auto dummy = thrust::make_tuple(fpcomplex(0.,0.), fpcomplex(0.,0.));
                
                integrals[i*n_res + j] = thrust::transform_reduce(
                    thrust::make_zip_iterator(thrust::make_tuple(binIndex, arrayAddress, effFunc)),
                    thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, arrayAddress, effFunc)),
                    *(integrators[i][j]),
                    dummy,
                    reduce_func);

                auto int_i = thrust::get<0>(integrals[i*n_res + j]);

                // if(i==j)
                //     std::cout << "Integral " << i << j << "= "<< int_i.real()*binSizeFactor << "," << int_i.imag()*binSizeFactor <<  "\n";
                // if(i<j)
                //     std::cout << "Integral " << i << j << "= "<< 2.*int_i.real()*binSizeFactor << "," << 2.*int_i.imag()*binSizeFactor <<  "\n";
         
        }
    }
    cudaDeviceSynchronize();

    for(int i = 0; i < n_res; ++i){
        
        fpcomplex int_i = thrust::get<0>(integrals[i*n_res + i])*binSizeFactor; //int |F_i|^2 dmdth
        calculators[i]->setResonanceIndex(decayInfo.resonances[i]->getFunctionIndex());
        calculators[i]->setDalitzIndex(getFunctionIndex());
        calculators[i]->setNorm(int_i.real());
        if(redoIntegral[i]){
            //std::cout << "Integral " << decayInfo.resonances[i]->getName() << "= "<< int_i.real() << "," << int_i.imag() <<  "\n";
            thrust::transform(
                        thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
                        thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
                        strided_range<thrust::device_vector<fpcomplex>::iterator>(cachedWaves[i]->begin(), cachedWaves[i]->end(), 1).begin(),
                        *(calculators[i])
            );
        }
    }

    cudaDeviceSynchronize();

    // End of time-consuming integrals.
    fpcomplex sumIntegral(0, 0);
    for(unsigned int i = 0; i < n_res; ++i) {
        // int param_i = parameters + resonanceOffset_DP + resonanceSize * i;
        fptype mag = host_parameters[parametersIdx + i * 2 + 1];
        fptype phs = host_parameters[parametersIdx + i * 2 + 2];

        const fpcomplex amplitude_i = thrust::polar(mag,phs);

        for(unsigned int j = 0; j < n_res; ++j) {
            // int param_j = parameters + resonanceOffset_DP + resonanceSize * j;
            mag = host_parameters[parametersIdx + j * 2 + 1];
            phs = -host_parameters[parametersIdx + j * 2 + 2];
            const fpcomplex amplitude_j = thrust::polar(mag,phs);

            fptype int_i = thrust::get<0>(integrals[i*n_res + i]).real();
            fptype int_j = thrust::get<0>(integrals[j*n_res + j]).real();

            if(int_i<1e-10 || int_j<1e-10)
                return 0;

            fptype fNorm_i = 1./int_i; // get<0> returns integral of RBW widthout eff
            fptype fNorm_j = 1./int_j;
            fpcomplex int_ij = thrust::get<1>(integrals[i*n_res + j]); // int_ij= integral_{ij}*eff

            sumIntegral += amplitude_i*amplitude_j*int_ij*sqrt(fNorm_i)*sqrt(fNorm_j);
        }
    }

    fptype ret           = sumIntegral.real(); // That complex number is a square, so it's fully real
    host_normalizations[normalIdx + 1] = 1.0 / ret;
    cachedNormalization                = 1.0 / ret;
    return ret;
}

__host__ auto Amp3Body::sumCachedWave(size_t i) const -> fpcomplex {
    const thrust::device_vector<fpcomplex> &vec = getCachedWaveNoCopy(i);

    fpcomplex ret = thrust::reduce(vec.begin(), vec.end(), fpcomplex(0, 0), thrust::plus<fpcomplex>());

    return ret;
}

__host__ auto Amp3Body::getCachedWave(size_t i) const -> const std::vector<std::complex<fptype>> {
    // TODO: This calls itself immediately ?
    auto ret_thrust = getCachedWave(i);
    std::vector<std::complex<fptype>> ret(ret_thrust.size());
    thrust::copy(ret_thrust.begin(), ret_thrust.end(), ret.begin());
    return ret;
}

__host__ auto Amp3Body::fit_fractions(bool print, std::string print_to_file_path) -> std::vector<std::vector<fptype>> {
   
    recursiveSetNormalization(1.0);

    host_normalizations.sync(d_normalizations);
    
    size_t n_res     = getDecayInfo().resonances.size();
    size_t totalBins = _m13.getNumBins() * _m23.getNumBins();
    double binSizeFactor = 1;
    binSizeFactor *= _m13.getBinSize();
    binSizeFactor *= _m23.getBinSize();

    if(!dalitzNormRange) {
        gooMalloc((void **)&dalitzNormRange, 6 * sizeof(fptype));
    }

    // This line runs once
    static std::array<fptype, 6> host_norms{{0, 0, 0, 0, 0, 0}};

    std::array<fptype, 6> current_host_norms{{_m13.getLowerLimit(),
                                              _m13.getUpperLimit(),
                                              static_cast<fptype>(_m13.getNumBins()),
                                              _m23.getLowerLimit(),
                                              _m23.getUpperLimit(),
                                              static_cast<fptype>(_m23.getNumBins())}};

    if(host_norms != current_host_norms) {
        host_norms = current_host_norms;
    }

    MEMCPY(dalitzNormRange, host_norms.data(), 6 * sizeof(fptype), cudaMemcpyHostToDevice);


    // Only do this bit if masses or widths have changed.
   thrust::constant_iterator<fptype *> arrayAddress(dalitzNormRange);
    thrust::counting_iterator<int> binIndex(0);
   
    thrust::constant_iterator<fptype *> dataArray(dev_event_array);
    thrust::constant_iterator<int> eventSize(totalEventSize);
    thrust::counting_iterator<int> eventIndex(eventOffset);


     auto reduce_func = [] __host__ __device__ (const thrust::tuple<fpcomplex, fpcomplex>& x, const thrust::tuple<fpcomplex, fpcomplex>& y)
    {
        return thrust::make_tuple(thrust::get<0>(x) + thrust::get<0>(y), thrust::get<1>(x) + thrust::get<1>(y));
    };

    for(int i = 0; i < n_res; ++i) {
        for(int j = 0; j < n_res; ++j) {

            integrators[i][j]->setDalitzIndex(getFunctionIndex());
            integrators[i][j]->setResonanceIndex(decayInfo.resonances[i]->getFunctionIndex());
            integrators[i][j]->setEfficiencyIndex(decayInfo.resonances[j]->getFunctionIndex());
            thrust::constant_iterator<int> effFunc(efficiencyFunction);
            auto dummy = thrust::make_tuple(fpcomplex(0.,0.), fpcomplex(0.,0.));
            
            integrals[i*n_res + j] = thrust::transform_reduce(
                   thrust::make_zip_iterator(thrust::make_tuple(binIndex, arrayAddress, effFunc)),
                thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, arrayAddress, effFunc)),
                *(integrators[i][j]),
                dummy,
                reduce_func);

            // printf("%d %d \n",i,j);
        }
    }

    cudaDeviceSynchronize();

    // // End of time-consuming integrals.
    fpcomplex sumIntegral(0, 0);
    std::vector<std::vector<fptype>> AmpFFs(n_res, std::vector<fptype>(n_res));
    std::vector<std::vector<fptype>> AmpIntegral(n_res, std::vector<fptype>(n_res));
    
    for(unsigned int i = 0; i < n_res; ++i) {
        fptype mag = host_parameters[parametersIdx + i * 2 + 1];
        fptype phs = host_parameters[parametersIdx + i * 2 + 2];
        const fpcomplex amplitude_i = thrust::polar(mag,phs);
        fpcomplex buffer(0.,0.);

        for(unsigned int j = 0; j < n_res; ++j) {
            mag = host_parameters[parametersIdx + j * 2 + 1];
            phs = -host_parameters[parametersIdx + j * 2 + 2];
            const fpcomplex amplitude_j = thrust::polar(mag,phs);

            //fpcomplex amp_prod = amplitude_i*amplitude_j;

            fptype fNorm_i = 1./(thrust::get<0>(integrals[i*n_res + i]).real());
            fptype fNorm_j = 1./(thrust::get<0>(integrals[j*n_res + j]).real());
            fpcomplex int_ij = thrust::get<0>(integrals[i*n_res + j]);

            // if(i==j)
            //      buffer = thrust::norm(amplitude_i*int_ij.real());
            // else
            //      buffer = 2.*(amp_prod*int_ij).real();

            buffer = amplitude_i*amplitude_j*int_ij*sqrt(fNorm_i)*sqrt(fNorm_j);


            AmpFFs[i][j] = buffer.real() ;
            AmpIntegral[i][j] = int_ij.real()*binSizeFactor;
            sumIntegral += buffer;
        }
    }
    std::cout << "\n";

    totalFF_integral = sumIntegral.real();

    std::cout << "totalFF_integral (%)" <<totalFF_integral << std::endl;

    for(int i = 0; i < n_res; i++) {
        for(int j = 0; j < n_res; j++) {
            AmpFFs[i][j] /= totalFF_integral;
            AmpFFs[i][j] *= 100;
            if(i<j) AmpFFs[i][j] *= 2.;
            if(i>j) AmpFFs[i][j] *= 0.;

        }
    }

  
   
    Eigen::MatrixXd m(n_res, n_res);
    Eigen::MatrixXd m_int(n_res, n_res);
    for(int i = 0; i < n_res; i++){
        m.row(i) = Eigen::Map<Eigen::VectorXd>(&AmpFFs[i][0], n_res);
        m_int.row(i) = Eigen::Map<Eigen::VectorXd>(&AmpIntegral[i][0], n_res);
    }
    if(print){
        std::cout << "Fit Fractions Matrix (%): \n";
        std::cout << "*Note: the order of diag FFs is equal to the order that which resonances are pushed into the "
                        "resonance vector. \n";
        std::cout << std::fixed << m << std::endl;
    }
  
    fptype sumdiagffs = 0.;

    if(print){
        std::cout << "\n ";
        std::cout << "Diagonal Fit Fractions (%): \n";


        for(int i = 0; i < n_res; i++){
            auto name = decayInfo.resonances[i]->getName();
            std::cout  << name << "\t" << std::fixed << m(i, i) << '\n';
            sumdiagffs += m(i, i);
        }
        std::cout << "Sum of Diag FFs: " << sumdiagffs << "\n";
        std::cout << "\n";
    }

    if(print){
        std::cout << "Res Integral Matrix: \n";
        std::cout << std::fixed << (SymDp ? 2.*m_int :  m_int) << std::endl;
    }

    if(!print_to_file_path.empty()){
        std::ofstream file(print_to_file_path);
        file << "Fit Fractions Matrix (%): \n";
        file << "*Note: the order of diag FFs is equal to the order that which resonances are pushed into the "
                "resonance vector. \n";
        file << std::fixed << m << std::endl;
        file << "\n ";
        file << "Diagonal Fit Fractions (%): \n";

        for(int i = 0; i < n_res; i++){
            auto name = decayInfo.resonances[i]->getName();
            file << name << "\t" << std::fixed << m(i, i) << '\n';
        }
        file << "Sum of Diag FFs: " << sumdiagffs << "\n";
        file << "\n";

        std::cout << "Res Integral Matrix: \n";
        file << std::fixed << (SymDp ? 2.*m_int :  m_int) << std::endl;
        file << "Total integral: " << std::fixed << totalFF_integral << std::endl;
        file.close();
    }

    return AmpFFs;
}

__host__ auto Amp3Body::GenerateSig(unsigned int numEvents, int seed) -> std::
    tuple<mcbooster::ParticlesSet_h, mcbooster::VariableSet_h, mcbooster::RealVector_h, mcbooster::RealVector_h> {
    // Must configure our functions before any calculations!
    // setupObservables();
    // setIndices();

    initialize();

    // Defining phase space
    std::vector<mcbooster::GReal_t> masses{decayInfo.daug1Mass, decayInfo.daug2Mass, decayInfo.daug3Mass};
    mcbooster::PhaseSpace phsp(decayInfo.motherMass, masses, numEvents, generation_offset);

    if(seed != 0) {
        phsp.SetSeed(seed);
    } else {
        GOOFIT_INFO("Current generator seed {}, offset {}", phsp.GetSeed(), generation_offset);
    }

    // Generating numEvents events. Events are all generated inside the phase space with uniform distribution in
    // momentum space. Events must be weighted to have phase space distribution
    phsp.Generate(mcbooster::Vector4R(decayInfo.motherMass, 0.0, 0.0, 0.0));

    phsp.Unweight();

    auto d1 = phsp.GetDaughters(0);
    auto d2 = phsp.GetDaughters(1);
    auto d3 = phsp.GetDaughters(2);

    mcbooster::ParticlesSet_d pset(3);
    pset[0] = &d1;
    pset[1] = &d2;
    pset[2] = &d3;

    auto SigGen_M12_d = mcbooster::RealVector_d(numEvents);
    auto SigGen_M13_d = mcbooster::RealVector_d(numEvents);
    auto SigGen_M23_d = mcbooster::RealVector_d(numEvents);

    mcbooster::VariableSet_d VarSet_d(3);
    VarSet_d[0] = &SigGen_M12_d;
    VarSet_d[1] = &SigGen_M23_d;
    VarSet_d[2] = &SigGen_M13_d;

    // Evaluating invariant masses for each event
    Dim2 eval = Dim2();
    mcbooster::EvaluateArray<Dim2>(eval, pset, VarSet_d);

    mcbooster::VariableSet_d GooVarSet_d(3);
    GooVarSet_d[0] = VarSet_d[0];
    GooVarSet_d[1] = VarSet_d[2];
    GooVarSet_d[2] = VarSet_d[1];

    auto h1 = new mcbooster::Particles_h(d1);
    auto h2 = new mcbooster::Particles_h(d2);
    auto h3 = new mcbooster::Particles_h(d3);

    mcbooster::ParticlesSet_h ParSet(3);
    ParSet[0] = h1;
    ParSet[1] = h2;
    ParSet[2] = h3;

    auto SigGen_M12_h = new mcbooster::RealVector_h(SigGen_M12_d);
    auto SigGen_M23_h = new mcbooster::RealVector_h(SigGen_M23_d);
    auto SigGen_M13_h = new mcbooster::RealVector_h(SigGen_M13_d);

    mcbooster::VariableSet_h VarSet(3);
    VarSet[0] = SigGen_M12_h;
    VarSet[1] = SigGen_M23_h;
    VarSet[2] = SigGen_M13_h;

    mcbooster::RealVector_d weights(phsp.GetWeights());
    phsp.FreeResources();

    auto DS = new mcbooster::RealVector_d(3 * numEvents);
    thrust::counting_iterator<int> eventNumber(0);

#pragma unroll

    for(int i = 0; i < 2; ++i) {
        mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr(DS->begin() + i, DS->end(), 3);
        thrust::copy(GooVarSet_d[i]->begin(), GooVarSet_d[i]->end(), sr.begin());
    }

    mcbooster::strided_range<mcbooster::RealVector_d::iterator> sr(DS->begin() + 2, DS->end(), 3);
    thrust::copy(eventNumber, eventNumber + numEvents, sr.begin());

    // Giving events to GooFit. Format of dev_evt_array must be (s12, s13, eventNumber). s23 is calculated automatically
    // in src/PDFs/physics/detail/SpecialResonanceCalculator.cu
    dev_event_array = thrust::raw_pointer_cast(DS->data());
    setDataSize(numEvents, 3);

    generation_no_norm = true; // we need no normalization for generation, but we do need to make sure that norm = 1;
    SigGenSetIndices();
    copyParams();
    normalize();
    setForceIntegrals();
    host_normalizations.sync(d_normalizations);

    auto fc = fitControl;
    setFitControl(std::make_shared<ProbFit>());

    thrust::device_vector<fptype> results;
    GooPdf::evaluate_with_metric(results);

    // evaluating amplitudes for generated events, amplitudes are incorporated in weights
    thrust::transform(
        results.begin(), results.end(), weights.begin(), weights.begin(), thrust::multiplies<mcbooster::GReal_t>());

    // Filing accept/reject flags for resonant distribution for each generated event
    mcbooster::BoolVector_d flags(numEvents);
    fillMCFlags(flags, weights, numEvents);

    auto weights_h = mcbooster::RealVector_h(weights);
    auto results_h = mcbooster::RealVector_h(results);
    auto flags_h   = mcbooster::BoolVector_h(flags);
    cudaDeviceSynchronize();

    setFitControl(fc);

    return std::make_tuple(ParSet, VarSet, weights_h, flags_h);
}
} // namespace GooFit
