#include <goofit/PDFs/basic/poissonPdf.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <math.h>

//using namespace std;

namespace GooFit {
	__device__ fptype device_PoissonPdf(fptype *evt, ParameterContainer &pc) {
		int id = pc.getObservable(0);
		fptype lambda = pc.getParameter(0);
		fptype x = RO_CACHE(evt[id]);
		pc.incrementIndex(1,1,0,1,1);
		fptype ret =  exp(-1*lambda)*pow(lambda,x)/tgamma(x+1); //doubt this will be final return statement, do I need to define another function to get this to work? I don't know if tgamma works, it might not.
		return ret;
	}


__device__ device_function_ptr ptr_to_PoissonPdf = device_PoissonPdf;

__host__ PoissonPdf::PoissonPdf(std::string n, Observable _x, Variable lambda) : GooPdf("PoissonPdf",n,_x,lambda) {
	registerFunction("ptr_to_PoissonPdf",ptr_to_PoissonPdf);
	initialize();
} 
}
