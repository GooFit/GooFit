#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/TrigThresholdPdf.h>

namespace GooFit {

__device__ auto threshCalc(fptype distance, fptype linConst) -> fptype {
    fptype ret
        = (distance > fptype(0.5) ? fptype(1) : (linConst + (1 - linConst) * sin(distance * fptype(3.14159265))));
    return ret;
}

__device__ auto device_TrigThresholdUpper(fptype *evt, ParameterContainer &pc) -> fptype {
    int id           = pc.getObservable(0);
    fptype x         = RO_CACHE(evt[id]);
    fptype thresh    = pc.getParameter(0);
    fptype trigConst = pc.getParameter(1);
    fptype linConst  = pc.getParameter(2);

    pc.incrementIndex(1, 3, 0, 1, 1);

    trigConst *= (thresh - x);
    return threshCalc(trigConst, linConst);
}

__device__ auto device_TrigThresholdLower(fptype *evt, ParameterContainer &pc) -> fptype {
    int id           = pc.getObservable(0);
    fptype x         = RO_CACHE(evt[id]);
    fptype thresh    = pc.getParameter(0);
    fptype trigConst = pc.getParameter(1);
    fptype linConst  = pc.getParameter(2);

    pc.incrementIndex(1, 3, 0, 1, 1);

    trigConst *= (x - thresh);
    return threshCalc(trigConst, linConst);
}

__device__ auto device_VerySpecialEpisodeTrigThresholdUpper(fptype *evt, ParameterContainer &pc) -> fptype {
    // Annoying special case for use with Mikhail's efficiency function across the Dalitz plot
    int id_x = pc.getObservable(0);
    int id_y = pc.getObservable(1);

    fptype x = RO_CACHE(evt[id_x]);
    fptype y = RO_CACHE(evt[id_y]);

    fptype thresh    = pc.getParameter(0);
    fptype trigConst = pc.getParameter(1);
    fptype linConst  = pc.getParameter(2);
    fptype z         = pc.getParameter(3) - x - y;

    pc.incrementIndex(1, 4, 0, 2, 1);

    trigConst *= (thresh - z);
    return threshCalc(trigConst, linConst);
}

__device__ auto device_VerySpecialEpisodeTrigThresholdLower(fptype *evt, ParameterContainer &pc) -> fptype {
    int id_x = pc.getObservable(0);
    int id_y = pc.getObservable(1);

    fptype x = RO_CACHE(evt[id_x]);
    fptype y = RO_CACHE(evt[id_y]);

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
    : GooPdf("TrigThresholdPdf", n, _x, thresh, trigConst, linConst) {
    if(upper)
        registerFunction("ptr_to_TrigThresholdUpper", ptr_to_TrigThresholdUpper);
    else
        registerFunction("ptr_to_TrigThresholdLower", ptr_to_TrigThresholdLower);

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
    : GooPdf("TrigThresholdPdf", n, _x, _y, thresh, trigConst, linConst, massConstant) {
    if(upper)
        registerFunction("ptr_to_VerySpecialEpisodeTrigThresholdUpper", ptr_to_VerySpecialEpisodeTrigThresholdUpper);
    else
        registerFunction("ptr_to_VerySpecialEpisodeTrigThresholdLower", ptr_to_VerySpecialEpisodeTrigThresholdLower);

    initialize();
}

} // namespace GooFit
