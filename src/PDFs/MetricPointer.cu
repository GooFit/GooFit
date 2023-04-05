#include <goofit/Error.h>
#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/MetricPointer.h>
#include <goofit/PDFs/detail/Globals.h>

namespace GooFit {

__device__ auto calculateEval(fptype rawPdf, fptype *evtVal, fptype norm) -> fptype {
    // Just return the raw PDF value, for use in (eg) normalization.
    return rawPdf;
}

__device__ auto calculateNLL(fptype rawPdf, fptype *evtVal, fptype norm) -> fptype {
    rawPdf *= norm;
    return rawPdf > 0.0 ? -log(rawPdf) : 0.0;
}

__device__ auto calculateProb(fptype rawPdf, fptype *evtVal, fptype norm) -> fptype {
    // Return probability, ie normalized PDF value.
    return rawPdf * norm;
}

__device__ auto calculateBinAvg(fptype rawPdf, fptype *evtVal, fptype norm) -> fptype {
    // TODO:(brad) address these metric devices later
    rawPdf *= norm;
    rawPdf *= evtVal[1]; // Bin volume

    // Log-likelihood of numEvents with expectation of exp is (-exp + numEvents*ln(exp) - ln(numEvents!)).
    // The last is constant, so we drop it; and then multiply by minus one to get the negative log-likelihood.
    if(rawPdf > 0) {
        fptype expEvents = c_totalEvents * rawPdf;
        return (expEvents - evtVal[0] * log(expEvents));
    }

    return 0;
}

__device__ auto calculateBinWithError(fptype rawPdf, fptype *evtVal, fptype norm) -> fptype {
    // TODO:(brad) address these metric devices later

    // In this case interpret the rawPdf as just a number, not a number of events.
    // Do not divide by integral over phase space, do not multiply by bin volume,
    // and do not collect 200 dollars. evtVal should have the structure (bin entry, bin error).
    // printf("[%i, %i] ((%f - %f) / %f)^2 = %f\n", BLOCKIDX, THREADIDX, rawPdf, evtVal[0], evtVal[1], pow((rawPdf -
    // evtVal[0]) / evtVal[1], 2));
    rawPdf -= evtVal[0]; // Subtract observed value.
    rawPdf /= evtVal[1]; // Divide by error.
    rawPdf *= rawPdf;
    return rawPdf;
}

__device__ auto calculateChisq(fptype rawPdf, fptype *evtVal, fptype norm) -> fptype {
    // TODO:(brad) address these metric devices later
    rawPdf *= norm;
    rawPdf *= evtVal[1]; // Bin volume

    return POW2(rawPdf * c_totalEvents - evtVal[0]) / (evtVal[0] > 1 ? evtVal[0] : 1);
}

__device__ device_metric_ptr ptr_to_Eval         = calculateEval;
__device__ device_metric_ptr ptr_to_NLL          = calculateNLL;
__device__ device_metric_ptr ptr_to_Prob         = calculateProb;
__device__ device_metric_ptr ptr_to_BinAvg       = calculateBinAvg;
__device__ device_metric_ptr ptr_to_BinWithError = calculateBinWithError;
__device__ device_metric_ptr ptr_to_Chisq        = calculateChisq;

// 211222 mds functionPtrToNameMap is declared in Globals.cpp
// it is a map of the device pointer values to the names of the methods;
// meant to be used for debugging/tracking logic; it is not used
// in any calculations.

// fill the map each time a host_fcn_ptr is created
// so we can access the name later

auto getMetricPointer(EvalFunc val) -> void * {
    if(val == EvalFunc::Eval) {
        host_fcn_ptr                       = get_device_symbol_address(ptr_to_Eval);
        functionPtrToNameMap[host_fcn_ptr] = "calculateEval";
    } else if(val == EvalFunc::NLL) {
        host_fcn_ptr                       = get_device_symbol_address(ptr_to_NLL);
        functionPtrToNameMap[host_fcn_ptr] = "calculateNLL";
    } else if(val == EvalFunc::Prob) {
        host_fcn_ptr                       = get_device_symbol_address(ptr_to_Prob);
        functionPtrToNameMap[host_fcn_ptr] = "calculateProb";
    } else if(val == EvalFunc::BinAvg) {
        host_fcn_ptr                       = get_device_symbol_address(ptr_to_BinAvg);
        functionPtrToNameMap[host_fcn_ptr] = "calculateBinAvg";
    } else if(val == EvalFunc::BinWithError) {
        host_fcn_ptr                       = get_device_symbol_address(ptr_to_BinWithError);
        functionPtrToNameMap[host_fcn_ptr] = "calculateBinWithError";
    } else if(val == EvalFunc::Chisq) {
        host_fcn_ptr                       = get_device_symbol_address(ptr_to_Chisq);
        functionPtrToNameMap[host_fcn_ptr] = "calculateChisq";
    } else {
        throw GeneralError("Non-existent metric pointer choice");
    }
    GOOFIT_TRACE("Selecting {} for the metric pointer", evalfunc_to_string(val));

    return host_fcn_ptr;
}

} // namespace GooFit
