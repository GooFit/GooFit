#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/SqDpHistogramPdf.h>
#include <goofit/Variable.h>

#include <TRandom3.h>
#include <random>

namespace GooFit {

__constant__ fptype *dev_base_sqdp_histograms[100]; // Multiple histograms for the case of multiple PDFs

unsigned int SqDpHistogramPdf::totalHistograms = 0;

__device__ auto device_EvalSqDPHistogram(fptype *evt, ParameterContainer &pc) -> fptype {
    // Structure is
    // nP smoothingIndex totalHistograms (limit1 step1 bins1) (limit2 step2 bins2) nO o1 o2
    // where limit and step are indices into functorConstants.

    int numCons          = pc.getNumConstants();
    int numObs           = pc.getNumObservables();
    int numParms         = pc.getNumParameters();
    int numVars          = pc.getConstant(0);
    int globalBinNumber  = 0;
    int previous         = 1;
    int myHistogramIndex = pc.getConstant(1); // 1 only used for smoothing

    auto motherMass = pc.getConstant(numCons-5);
    auto d1Mass = pc.getConstant(numCons-4);
    auto d2Mass = pc.getConstant(numCons-3);
    auto d3Mass = pc.getConstant(numCons-2);
    bool isEffHist = pc.getConstant(numCons-1);
   

    int mp_Index = pc.getObservable(0);
    int thp_Index = pc.getObservable(1);

    auto mp = evt[mp_Index];
    auto thp = evt[thp_Index];

    auto jac = calc_SqDp_Jacobian(mp , thp, motherMass, d1Mass, d2Mass, d3Mass);

    for(int i = 0; i < numVars; ++i) {
        int varIndex      = pc.getObservable(i);
        int lowerBoundIdx = 3 * (i + 1);

        // don't use RO_CACHE as this is used as efficiency for Amp3Body
        fptype currVariable = evt[varIndex];
        fptype lowerBound   = pc.getConstant(i * 3 + 4);
        fptype step         = pc.getConstant(i * 3 + 5);

        currVariable -= lowerBound;
        currVariable /= step;

        auto localBinNumber = static_cast<int>(floor(currVariable));
        globalBinNumber += previous * localBinNumber;

        // This is accessing too far ahead?
        int offset = pc.getConstant(lowerBoundIdx);
        previous *= offset;
    }
    // printf("histIndex = %d \n",myHistogramIndex);
    
    fptype *myHistogram = dev_base_sqdp_histograms[myHistogramIndex];
    fptype ret          = myHistogram[globalBinNumber];

    pc.incrementIndex(1, numParms, numCons, numObs, 1);

    // if(isEffHist)
    //     jac = 1.;
  

    // return ret;

    // if(isEffHist)
    //     return ret*=jac;
    // else 
    //     return ret;
    
    const fptype D0MASS       = 1.865;
    const fptype range        = 0.030;
    const fptype D0veto_min   = (D0MASS - range) * (D0MASS - range);
    const fptype D0veto_max   = (D0MASS + range) * (D0MASS + range);

    
    auto m12 = calc_m12(mp,motherMass, d1Mass, d2Mass, d3Mass);
    auto m13 = calc_m13(m12,cos(M_PI*thp),motherMass, d1Mass, d2Mass, d3Mass);
    auto s23 = motherMass*motherMass + d1Mass*d1Mass +  d2Mass*d2Mass + d3Mass*d3Mass - m12*m12 - m13*m13;
    auto s13 = m13*m13;

    if((s13<D0veto_min || s13>D0veto_max) && (s23<D0veto_min || s23>D0veto_max)){
        if(isEffHist)
            return ret*=jac;
        else 
            return ret;
    }else{
        return 0.0;
    }
   
}

__device__ device_function_ptr ptr_to_EvalSqDPHistogram = device_EvalSqDPHistogram;

__host__ SqDpHistogramPdf::SqDpHistogramPdf(std::string n, TH2* hist, Observable &mPrime, Observable &thPrime, Variable motherMass, Variable d1Mass, Variable d2Mass, Variable d3Mass, bool upperHalfOnly, bool eff)
    : GooPdf("SqDpHistogramPdf", n, mPrime, thPrime, motherMass, d1Mass, d2Mass, d3Mass) {
    
    mPrime.setNumBins(hist->GetNbinsX());
    thPrime.setNumBins(hist->GetNbinsY());

    BinnedDataSet *binDataSet     = new BinnedDataSet({mPrime,thPrime});

    fptype dx = abs(mPrime.getLowerLimit()-mPrime.getUpperLimit())/mPrime.getNumBins();
    fptype dy = abs(thPrime.getLowerLimit()-thPrime.getUpperLimit())/thPrime.getNumBins();

    // //fill binDataSet
    for(int i=0; i<mPrime.getNumBins();i++){
        for(int j=0; j<thPrime.getNumBins();j++){
           
            fptype mp = mPrime.getLowerLimit() + (i+0.5)*(mPrime.getUpperLimit() - mPrime.getLowerLimit())/(mPrime.getNumBins());
            fptype thp = thPrime.getLowerLimit() + (j+0.5)*(thPrime.getUpperLimit() - thPrime.getLowerLimit())/(thPrime.getNumBins());

            mPrime.setValue(mp);
            thPrime.setValue(thp);
    
            if (upperHalfOnly){
                if(thp>0.5)
                    thp = 1. - thp;
            }

            auto currBin = hist->FindBin(mp,thp);
            // fptype currVal = hist->GetBinContent(currBin);
            fptype currVal = hist->Interpolate(mp,thp);
            if(currVal<0.)
                continue;
                
            // if(eff)
            //     currVal = currVal>1. ? 1.: currVal; //protect > 1 eff
               
            binDataSet->addWeightedEvent(currVal);
        }                
    }

    std::cout << "BinDataSet filled with " << binDataSet->getNumEvents() << std::endl;

    int numVars = binDataSet->numVariables();
    int varIndex = 0;

    registerConstant(binDataSet->getObservables().size());
    registerConstant(totalHistograms);
    registerConstant(numVars);

    for(const Observable &var : binDataSet->getObservables()) {
        registerObservable(var);
        registerConstant(var.getNumBins());
        registerConstant(var.getLowerLimit());
        registerConstant(var.getBinSize());
        varIndex++;
    }

    registerConstant(motherMass.getValue());
    registerConstant(d1Mass.getValue());
    registerConstant(d2Mass.getValue());
    registerConstant(d3Mass.getValue());
    registerConstant(eff);
    

    unsigned int numbins = binDataSet->getNumBins();
    thrust::host_vector<fptype> host_histogram;

    fptype integral = 0.;

    for(unsigned int i = 0; i < numbins; ++i) {
        fptype curr = binDataSet->getBinContent(i);
        host_histogram.push_back(curr);
        integral+=curr;
    }

    std::cout << "host_histogram filled with " << host_histogram.size() << std::endl;

    integral*=dx*dy;
   
    // printf("Integral*dx*dy %f \n",totalEvents);
//    if(!eff){
//         thrust::transform(host_histogram.begin(),
//                     host_histogram.end(), 
//                     thrust::make_constant_iterator(1./totalEvents),
//                     host_histogram.begin(),
//                     thrust::multiplies<fptype>());
//     }

    std::cout << "total " << getName() << " " << integral << std::endl;

    if(integral > 0){
        copyHistogramToDevice(host_histogram);
    }else{
        std::cout << "Warning: Empty histogram supplied to " << getName()
                  << " not copied to device. Expect copyHistogramToDevice call later.\n";
    }
    registerFunction("ptr_to_EvalSqDPHistogram", ptr_to_EvalSqDPHistogram);

    initialize();
}

auto pointerToFirstHist(thrust::device_vector<fptype> *hist) -> fptype * { return (&((*hist)[0])).get(); }

auto pointerToFirstHist(thrust::host_vector<fptype> *hist) -> fptype * {
    // (*hist) is the host_vector.
    // (*hist)[0] is a 'reference' - Thrust class, not ordinary C++ reference -
    // to the first element of the vector.
    // &((*hist)[0]) is a 'Pointer', as defined by the host_vector, to the location
    // of the 'reference'. Fortunately this is by default fptype*!
    return &((*hist)[0]);
}

__host__ void SqDpHistogramPdf::copyHistogramToDevice(thrust::host_vector<fptype> &host_histogram) {

    dev_base_sqdp_histogram     = new thrust::device_vector<fptype>(host_histogram);
  
    static fptype *dev_address[1] ;
    dev_address[0]= pointerToFirstHist(dev_base_sqdp_histogram);

    MEMCPY_TO_SYMBOL(
        dev_base_sqdp_histograms, dev_address, sizeof(fptype *), totalHistograms * sizeof(fptype *), cudaMemcpyHostToDevice);
  

    totalHistograms++;

     std::cout << "totalHistograms = " << totalHistograms << std::endl;

    int expectedBins = 1;

    for(Observable &observable : observablesList) {
        expectedBins *= observable.getNumBins();
    }

    if(expectedBins != host_histogram.size()) {
        std::cout << "Warning: Host Histogram supplied to " << getName() << " has " << host_histogram.size()
                  << " bins, expected " << expectedBins << " - may indicate a problem.\n";
    }

  
}


__host__ auto SqDpHistogramPdf::normalize() -> fptype {

    //return totalEvents;
  
    fptype ret = thrust::reduce(dev_base_sqdp_histogram->begin(), dev_base_sqdp_histogram->end()); 

     for(unsigned int varIndex = 0; varIndex < observablesList.size(); ++varIndex) {
        fptype binSize = constantsList[3 + 3 * varIndex + 2];
        ret *= binSize; // Bin size cached by constructor.
    }
   
   

    host_normalizations[normalIdx + 1] = 1.0 / ret;
    cachedNormalization                = 1.0 / ret;

    return ret;
}
} // namespace GooFit
