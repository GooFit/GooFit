#include <goofit/PDFs/basic/chebyshevPdf.h>
#include <goofit/PDFs/ParameterContainer.h>


namespace GooFit {

	__device__ fptype chebyshev(int max, fptype x) {

	fptype s = 0.0;
	switch(max) { //My idea: use a switch case to determine whether to call itself recursively or degrees 0 / 1 values
		case 0: //T(0,x) = 1; Hard coded expression
			s += 1.0;
			break;
		case 1: //T(1,x) = x; Other necessary hard-coded expression
			s+= x;
			break;
//		case (n < 0): //To catch people who input negative integers.
//			cout << "Not a valid value for n, try again: ";
//			break;
		default: //T(n,x) = 2xT(n-1,x) - T(n-2,x) for all max >= 2.
			s += (2*x*chebyshev(max-1,x))-chebyshev(max-2,x);
	}
	return s;
}


	__device__ fptype device_chebyshevPdf(fptype *evt, ParameterContainer &pc) {
		int id       = pc.getObservable(0);
		fptype ret = 0;
		int max = pc.getConstant(0);         //What'll the max variable be? Do I get it from pc?
		fptype x     = RO_CACHE(evt[id]); //The legendre fitting, I think, starts from PolynomialPdf
		//From there, it goes backwards to determine the closest sum of chebyshev Polynomials.
		/*But the Polynomial fit uses an array of coefficients, how will that be passed? 
		  Do I need to change the input parameters? I'll ask Daniel
		 */
		for(int j = 0; j <= max; j++) { //What'll the loop limit be? 
			fptype p = pc.getParameter(j);
			ret += p*chebyshev(x,j);
		}
		pc.incrementIndex(1, max, 1, 1, 1);
		return ret;
	}

	__device__ device_function_ptr ptr_to_chebyshevPdf = device_chebyshevPdf;

	__host__ chebyshevPdf::chebyshevPdf(std::string n, Observable _x, std::vector<Variable> weights,  unsigned int max)
		: GooPdf("chebyshevPdf", n, _x) { //Do I need 6 or 5 params?
			registerFunction("ptr_to_chebyshev", ptr_to_chebyshevPdf);
			registerConstant(max);
			for(Variable &v : weights) {
				registerParameter(v);
			}
			initialize();
		}

}
