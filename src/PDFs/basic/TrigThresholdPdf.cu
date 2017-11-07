#include "goofit/PDFs/basic/TrigThresholdPdf.h"

namespace GooFit {

__device__ fptype threshCalc(fptype distance, fptype linConst) {
    fptype ret
        = (distance > fptype(0.5) ? fptype(1) : (linConst + (1 - linConst) * sin(distance * fptype(3.14159265))));
    return ret;
}

__device__ fptype device_TrigThresholdUpper(fptype *evt, ParameterContainer &pc) {
    int id           = RO_CACHE(pc.observables[pc.observableIdx + 1]);
    fptype x         = evt[id];
    fptype thresh    = RO_CACHE(pc.parameters[pc.parameterIdx + 1]);
    fptype trigConst = RO_CACHE(pc.parameters[pc.parameterIdx + 2]);
    fptype linConst  = RO_CACHE(pc.parameters[pc.parameterIdx + 3]);

    pc.incrementIndex(1, 3, 0, 1, 1);

    trigConst *= (thresh - x);
    return threshCalc(trigConst, linConst);
}

__device__ fptype device_TrigThresholdLower(fptype *evt, ParameterContainer &pc) {
    int id           = RO_CACHE(pc.observables[pc.observableIdx + 2]);
    fptype x         = evt[id];
    fptype thresh    = RO_CACHE(pc.parameters[pc.parameterIdx + 1]);
    fptype trigConst = RO_CACHE(pc.parameters[pc.parameterIdx + 2]);
    fptype linConst  = RO_CACHE(pc.parameters[pc.parameterIdx + 3]);

    pc.incrementIndex(1, 3, 0, 1, 1);

    trigConst *= (x - thresh);
    return threshCalc(trigConst, linConst);
}

__device__ fptype device_VerySpecialEpisodeTrigThresholdUpper(fptype *evt, ParameterContainer &pc) {
    // Annoying special case for use with Mikhail's efficiency function across the Dalitz plot
    int id_x = RO_CACHE(pc.observables[pc.observableIdx + 1]);
    int id_y = RO_CACHE(pc.observables[pc.observableIdx + 2]);

    fptype x = evt[id_x];
    fptype y = evt[id_y];

    fptype thresh    = RO_CACHE(pc.parameters[pc.parameterIdx + 1]);
    fptype trigConst = RO_CACHE(pc.parameters[pc.parameterIdx + 2]);
    fptype linConst  = RO_CACHE(pc.parameters[pc.parameterIdx + 3]);
    fptype z         = RO_CACHE(pc.parameters[pc.parameterIdx + 4]) - x - y;

    pc.incrementIndex(1, 4, 0, 2, 1);

    trigConst *= (thresh - z);
    return threshCalc(trigConst, linConst);
}

__device__ fptype device_VerySpecialEpisodeTrigThresholdLower(fptype *evt, ParameterContainer &pc) {
    int id_x = RO_CACHE(pc.observables[pc.observableIdx + 1]);
    int id_y = RO_CACHE(pc.observables[pc.observableIdx + 2]);

    fptype x = evt[id_x];
    fptype y = evt[id_y];

    fptype thresh    = RO_CACHE(pc.parameters[pc.parameterIdx + 1]);
    fptype trigConst = RO_CACHE(pc.parameters[pc.parameterIdx + 2]);
    fptype linConst  = RO_CACHE(pc.parameters[pc.parameterIdx + 3]);
    fptype z         = RO_CACHE(pc.parameters[pc.parameterIdx + 4]) - x - y;

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
    std::vector<unsigned int> pindices;
    pindices.push_back(registerParameter(thresh));
    pindices.push_back(registerParameter(trigConst));
    pindices.push_back(registerParameter(linConst));

    if(upper) {
        GET_FUNCTION_ADDR(ptr_to_TrigThresholdUpper);
        trigThreshType = 0;
    } else {
        GET_FUNCTION_ADDR(ptr_to_TrigThresholdLower);
        trigThreshType = 1;
    }

    initialize(pindices);
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
    std::vector<unsigned int> pindices;
    pindices.push_back(registerParameter(thresh));
    pindices.push_back(registerParameter(trigConst));
    pindices.push_back(registerParameter(linConst));
    pindices.push_back(registerParameter(massConstant));

    if(upper) {
        GET_FUNCTION_ADDR(ptr_to_VerySpecialEpisodeTrigThresholdUpper);
        trigThreshType = 2;
    } else {
        GET_FUNCTION_ADDR(ptr_to_VerySpecialEpisodeTrigThresholdLower);
        trigThreshType = 3;
    }

    initialize(pindices);
}

void TrigThresholdPdf::recursiveSetIndices() {
    if(trigThreshType == 0) {
        GET_FUNCTION_ADDR(ptr_to_TrigThresholdUpper);
        GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_TrigThresholdUpper");
    } else if(trigThreshType == 1) {
        GET_FUNCTION_ADDR(ptr_to_TrigThresholdLower);
        GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName(), "ptr_to_TrigThresholdLower");
    } else if(trigThreshType == 2) {
        GET_FUNCTION_ADDR(ptr_to_VerySpecialEpisodeTrigThresholdUpper);
        GOOFIT_TRACE("host_function_table[{}] = {}({})",
                     num_device_functions,
                     getName(),
                     "ptr_to_VerySpecialEpisodeTrigThresholdUpper");
    } else if(trigThreshType == 3) {
        GET_FUNCTION_ADDR(ptr_to_VerySpecialEpisodeTrigThresholdLower);
        GOOFIT_TRACE("host_function_table[{}] = {}({})",
                     num_device_functions,
                     getName(),
                     "ptr_to_VerySpecialEpisodeTrigThresholdLower");
    }

    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx                               = num_device_functions++;

    populateArrays();
}

} // namespace GooFit
