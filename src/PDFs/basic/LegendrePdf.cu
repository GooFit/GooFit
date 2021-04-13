#include <goofit/PDFs/basic/Legendre.h>
#include <goofit/PDFs/ParameterContainer.h>


namespace GooFit {

	__device__ fptype Legendre(int max, fptype x) {

		int k = 0;
		fptype c = 1.0;
		fptype s = 0.0;
		c = 1/((pow(2,max))*Factorial(max));                    
		const int index = (2*max)+1;
		fptype p[index];

		for (int i=0;i <= max;i++)
		{

			p[2*i] = c*BinomCoeff(max,i)*(pow(-1,max-i));                           
			p[1+(2*i)] = 0;
		} //Expands polynomial out before taking the derivatives                         
		for (int j=max;j <= 2*max + 1; j++) 
		{
			p[j-max] = p[j]*(Factorial(j)/Factorial(j-max));
			if (max != 0) {
				p[j] = 0;
			} 	else {
				p[j] = 1; 
			}
		} //Calculates the coefficients to each polynomial term after all derivations
		while (k <= max) 
		{
			s += p[k]*(pow(x,k));
			//cout << p[k] << " ";
			k += 1;
		}
		return s;
	}

	__device__ int Factorial(int n) {
		if (n==0||n==1) {
			return 1;
		} 	else {
			return n*Factorial(n-1);
		}
	}

	__device__ fptype BinomCoeff(int n, int k) {
		return (Factorial(n)/(Factorial(k)*Factorial(n-k)));
	}

	__device__ fptype device_LegendrePdf(fptype *evt, ParameterContainer &pc) {
		int id       = pc.getObservable(0);
		fptype ret = 0;
		int max = pc.getConstant(0);         //What'll the max variable be? Do I get it from pc?
		printf("Constant is %d \n",pc.getConstant(0));
		printf("Offset is %d \n",pc.getParameter(0));
		fptype x     = RO_CACHE(evt[id]); //The legendre fitting, I think, starts from PolynomialPdf
		//From there, it goes backwards to determine the closest sum of Legendre Polynomials.
		/*But the Polynomial fit uses an array of coefficients, how will that be passed? 
		  Do I need to change the input parameters? I'll ask Daniel
		 */
		for(int j = 0; j <= max; j++) { //What'll the loop limit be? 
			fptype p = pc.getParameter(j);
			ret += p*Legendre(x,j);
		}
		pc.incrementIndex(1, max, 1, 1, 1);
		return ret;
	}

	__device__ device_function_ptr ptr_to_LegendrePdf = device_LegendrePdf;

	__host__ LegendrePdf::LegendrePdf(std::string n, Observable _x, std::vector<Variable> weights,  unsigned int max)
		: GooPdf("LegendrePdf", n, _x) { //Do I need 6 or 5 params?
			registerFunction("ptr_to_Legendre", ptr_to_LegendrePdf);
			registerConstant(max);
			for(Variable &v : weights) {
				registerParameter(v);
			}
			initialize();
		}

} // namespace GooFit
