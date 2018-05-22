#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/TrigThresholdPdf.h>

namespace GooFit {

__device__ fptype threshCalc(fptype distance, fptype linConst) {
    fptype ret
        = (distance > fptype(0.5) ? fptype(1) : (linConst + (1 - linConst) * sin(distance * fptype(3.14159265))));
    return ret;
}

__device__ fptype device_TrigThresholdUpper(fptype *evt, ParameterContainer &pc) {
    int id           = pc.getObservable(0);
    fptype x         = evt[id];
    fptype thresh    = pc.getParameter(0);
    fptype trigConst = pc.getParameter(1);
    fptype linConst  = pc.getParameter(2);

    pc.incrementIndex(1, 3, 0, 1, 1);

    trigConst *= (thresh - x);
    return threshCalc(trigConst, linConst);
}

__device__ fptype device_TrigThresholdLower(fptype *evt, ParameterContainer &pc) {
    int id           = pc.getObservable(0);
    fptype x         = evt[id];
    fptype thresh    = pc.getParameter(0);
    fptype trigConst = pc.getParameter(1);
    fptype linConst  = pc.getParameter(2);

    pc.incrementIndex(1, 3, 0, 1, 1);

    trigConst *= (x - thresh);
    return threshCalc(trigConst, linConst);
}

__device__ fptype device_VerySpecialEpisodeTrigThresholdUpper(fptype *evt, ParameterContainer &pc) {
    // Annoying special case for use with Mikhail's efficiency function across the Dalitz plot
    int id_x = pc.getObservable(0);
    int id_y = pc.getObservable(1);

    fptype x = evt[id_x];
    fptype y = evt[id_y];

    fptype thresh    = pc.getParameter(0);
    fptype trigConst = pc.getParameter(1);
    fptype linConst  = pc.getParameter(2);
    fptype z         = pc.getParameter(3) - x - y;

    pc.incrementIndex(1, 4, 0, 2, 1);

    trigConst *= (thresh - z);
    return threshCalc(trigConst, linConst);
}

__device__ fptype device_VerySpecialEpisodeTrigThresholdLower(fptype *evt, ParameterContainer &pc) {
    int id_x = pc.getObservable(0);
    int id_y = pc.getObservable(1);

    fptype x = evt[id_x];
    fptype y = evt[id_y];

    fptype thresh    = pc.getParameter(0);
    fptype trigConst = pc.getParameter(1);
    fptype linConst  = pc.getParameter(2);
    fptype z         = pc.getParameter(3) - x - y;

    pc.incrementIndex(1, 4, 0, 2, 1);

    trigConst *= (z - thresh);
    fptype ret = threshCalc(trigConst, linConst);

    // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
    // printf("TrigThreshold: (%f - %f = %f) -> %f %f\n", z, thresh, trigConst, linConst, ret);

    return ret;
}

__device__ device_function_ptr ptr_to_TrigThresholdUpper = device_TrigThresholdUpper;
__device__ device_function_ptr ptr_to_TrigThresholdLower = device_TrigThresholdLower;
__device__ device_function_ptr ptr_to_VerySpecialEpisodeTrigThresholdUpper
    = device_VerySpecialEpisodeTrigThresholdUpper;
__device__ device_function_ptr ptr_to_VerySpecialEpisodeTrigThresholdLower
    = device_VerySpecialEpisodeTrigThresholdLower;

__host__ TrigThresholdPdf::TrigThresholdPdf(
    std::string n, Observable _x, Variable thresh, Variable trigConst, Variable linConst, bool upper)
    : GooPdf(n, _x) {
    registerParameter(thresh);
    registerParameter(trigConst);
    registerParameter(linConst);

    if(upper) {
        trigThreshType = 0;
    } else {
        trigThreshType = 1;
    }

    initialize();
}

__host__ TrigThresholdPdf::TrigThresholdPdf(std::string n,
                                            Observable _x,
                                            Observable _y,
                                            Variable thresh,
                                            Variable trigConst,
                                            Variable linConst,
                                            Variable massConstant,
                                            bool upper)
    : GooPdf(n, _x, _y) {
    registerParameter(thresh);
    registerParameter(trigConst);
    registerParameter(linConst);
    registerParameter(massConstant);

    if(upper) {
        trigThreshType = 2;
    } else {
        trigThreshType = 3;
    }

    initialize();
}

void TrigThresholdPdf::recursiveSetIndices() {
    if(trigThreshType == 0) {
        GOOFIT_RECURSIVE_SET_INDICIES(ptr_to_TrigThresholdUpper);
    } else if(trigThreshType == 1) {
        GOOFIT_RECURSIVE_SET_INDICIES(ptr_to_TrigThresholdLower);
    } else if(trigThreshType == 2) {
        GOOFIT_RECURSIVE_SET_INDICIES(ptr_to_VerySpecialEpisodeTrigThresholdUpper);
    } else if(trigThreshType == 3) {
        GOOFIT_RECURSIVE_SET_INDICIES(ptr_to_VerySpecialEpisodeTrigThresholdLower);
    }
}

} // namespace GooFit
