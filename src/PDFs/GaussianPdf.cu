#include "goofit/PDFs/GaussianPdf.h"

EXEC_TARGET fptype device_Gaussian(fptype* evt, fptype* p, unsigned int* indices) {
    fptype x = evt[RO_CACHE(indices[2 + RO_CACHE(indices[0])])];
    fptype mean = RO_CACHE(p[RO_CACHE(indices[1])]);
    fptype sigma = RO_CACHE(p[RO_CACHE(indices[2])]);

    fptype ret = EXP(-0.5*(x-mean)*(x-mean)/(sigma*sigma));

    //if ((0 == THREADIDX) && (0 == BLOCKIDX)) cuPrintf("Gaussian Values %f %i %i %f %f %i\n", x, indices[1], indices[2], mean, sigma, callnumber);
    //cuPrintf("device_Gaussian %f %i %i %f %f %i %p %f\n", x, indices[1], indices[2], mean, sigma, callnumber, indices, ret);
    //if ((0 == THREADIDX) && (0 == BLOCKIDX))
    //printf("device_Gaussian %f %f %f %i %f\n", x, mean, sigma, callnumber, ret);


    return ret;
}

MEM_DEVICE device_function_ptr ptr_to_Gaussian = device_Gaussian;

__host__ GaussianPdf::GaussianPdf(std::string n, Variable* _x, Variable* mean, Variable* sigma)
    : GooPdf(_x, n) {
    std::vector<unsigned int> pindices;
    pindices.push_back(registerParameter(mean));
    pindices.push_back(registerParameter(sigma));
    GET_FUNCTION_ADDR(ptr_to_Gaussian);
    initialise(pindices);
}

__host__ fptype GaussianPdf::integrate(fptype lo, fptype hi) const {
    //static const fptype root2 = sqrt(2.);
    static const fptype rootPi = sqrt(atan2(0.0, -1.0));
    //static const fptype rootPiBy2 = rootPi / root2;

    unsigned int* indices = host_indices+parameters;
    //fptype xscale = root2*host_params[indices[2]];

    /*
    std::cout << "Gaussian integral: "
        << xscale << " "
        << host_params[indices[1]] << " "
        << host_params[indices[2]] << " "
        << ERF((hi-host_params[indices[1]])/xscale) << " "
        << ERF((lo-host_params[indices[1]])/xscale) << " "
        << rootPiBy2*host_params[indices[2]]*(ERF((hi-host_params[indices[1]])/xscale) -
    					  ERF((lo-host_params[indices[1]])/xscale))
        << std::endl;
    */
    //return rootPiBy2*host_params[indices[2]]*(ERF((hi-host_params[indices[1]])/xscale) -
    //					    ERF((lo-host_params[indices[1]])/xscale));

    // Integral over all R.
    fptype sigma = host_params[indices[2]];
    sigma *= root2*rootPi;
    return sigma;
}

