#include <goofit/Error.h>
#include <goofit/GlobalCudaDefines.h>
#include <goofit/PDFs/MetricPointer.h>
#include <goofit/PDFs/detail/Globals.h>

namespace GooFit {

__device__ fptype calculateEval(fptype rawPdf, fptype *evtVal, fptype norm) {
    // Just return the raw PDF value, for use in (eg) normalization.
    return rawPdf;
}

__device__ fptype calculateNLL(fptype rawPdf, fptype *evtVal, fptype norm) {
    rawPdf *= norm;
    return rawPdf > 0.0 ? -log(rawPdf) : 0.0;
}

__device__ fptype calculateProb(fptype rawPdf, fptype *evtVal, fptype norm) {
    // Return probability, ie normalized PDF value.
    return rawPdf * norm;
}

__device__ fptype calculateBinAvg(fptype rawPdf, fptype *evtVal, fptype norm) {
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

__device__ fptype calculateBinWithError(fptype rawPdf, fptype *evtVal, fptype norm) {
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

__device__ fptype calculateChisq(fptype rawPdf, fptype *evtVal, fptype norm) {
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

void *getMetricPointer(EvalFunc val) {
    if(val == EvalFunc::Eval) {
        host_fcn_ptr = get_device_symbol_address(ptr_to_Eval);
    } else if(val == EvalFunc::NLL) {
        host_fcn_ptr = get_device_symbol_address(ptr_to_NLL);
    } else if(val == EvalFunc::Prob) {
        host_fcn_ptr = get_device_symbol_address(ptr_to_Prob);
    } else if(val == EvalFunc::BinAvg) {
        host_fcn_ptr = get_device_symbol_address(ptr_to_BinAvg);
    } else if(val == EvalFunc::BinWithError) {
        host_fcn_ptr = get_device_symbol_address(ptr_to_BinWithError);
    } else if(val == EvalFunc::Chisq) {
        host_fcn_ptr = get_device_symbol_address(ptr_to_Chisq);
    } else {
        throw GeneralError("Non-existent metric pointer choice");
    }
    GOOFIT_TRACE("Selecting {} for the metric pointer", evalfunc_to_string(val));

    return host_fcn_ptr;
}

} // namespace GooFit
