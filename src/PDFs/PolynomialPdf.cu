#include "goofit/PDFs/basic/PolynomialPdf.h"
#include "goofit/Variable.h"
#include "goofit/Log.h"

namespace GooFit {

__device__ fptype device_Polynomial(fptype* evt, ParameterContainer &pc) {
    int id = RO_CACHE(pc.observables[pc.observableIdx + 1]);
    // Structure is nP lowestdegree c1 c2 c3 nO o1

    int numParams = RO_CACHE(pc.parameters[pc.parameterIdx]);
    int lowestDegree = RO_CACHE(pc.constants[pc.constantIdx + 1]);

    fptype x = evt[id];
    fptype ret = 0;

    //unsure why this starts at i=2...
    for(int i = 0; i < numParams; ++i) {
        fptype param = RO_CACHE(pc.parameters[pc.parameterIdx + i + 1]);
        ret += param * pow(x, lowestDegree + i);
    }

    pc.incrementIndex (1, numParams, 3, 1, 1);

    return ret;
}

__device__ fptype device_OffsetPolynomial(fptype* evt, ParameterContainer &pc) {
    int id = RO_CACHE(pc.observables[pc.observableIdx + 1]);

    int numParams = RO_CACHE(pc.parameters[pc.parameterIdx]);
    int lowestDegree = RO_CACHE(pc.constants[pc.constantIdx + 1]);

    fptype x = evt[id];
    //TODO: Not sure where this is pointing...
    //x -= RO_CACHE(p[RO_CACHE(indices[numParams])]);
    fptype ret = 0;

    for(int i = 2; i < numParams; ++i) {
        ret += RO_CACHE(pc.parameters[pc.parameterIdx + i]) * pow(x, lowestDegree + i - 2);
    }

    pc.incrementIndex(1, numParams, 1, 1, 1);

    return ret;
}

__device__ fptype device_MultiPolynomial(fptype* evt, ParameterContainer &pc) {
    int num_constants = RO_CACHE(pc.constants[pc.constantIdx]);
    int num_parameters = RO_CACHE(pc.parameters[pc.parameterIdx]);

    // Structure is nP, maxDegree, offset1, offset2, ..., coeff1, coeff2, ..., nO, o1, o2, ...

    int num_observables = RO_CACHE(pc.observables[pc.observableIdx]);
    int maxDegree = RO_CACHE(pc.constants[pc.constantIdx + 1]) + 1;
    // Only appears in construction (maxDegree + 1) or (x > maxDegree), so
    // may as well add the one and use >= instead.

    // Technique is to iterate over the full n-dimensional box, skipping matrix elements
    // whose sum of indices is greater than maxDegree. Notice that this is increasingly
    // inefficient as n grows, since a larger proportion of boxes will be skipped.
    int numBoxes = 1;

    for(int i = 0; i < num_observables; ++i)
        numBoxes *= maxDegree;

    int coeffNumber = num_observables; // Index of first coefficient is 2 + nO, not 1 + nO, due to maxDegree. (nO comes from offsets.)
    fptype ret = RO_CACHE(pc.parameters[pc.parameterIdx + 1 + coeffNumber]); // Coefficient of constant term.
    coeffNumber++;

    for(int i = 1; i < numBoxes;
        ++i) { // Notice skip of inmost 'box' in the pyramid, corresponding to all powers zero, already accounted for.
        fptype currTerm  = 1;
        int currIndex    = i;
        int sumOfIndices = 0;

        // if ((gpuDebug & 1) && (THREADIDX == 50) && (BLOCKIDX == 3))
        // if ((BLOCKIDX == internalDebug1) && (THREADIDX == internalDebug2))
        // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
        // printf("[%i, %i] Start box %i %f %f:\n", BLOCKIDX, THREADIDX, i, ret, evt[8]);
        for(int j = 0; j < num_observables; ++j) {
	    //TODO:Need to debug these
            int id = RO_CACHE(pc.observables[pc.observableIdx + 1 + j]);
            fptype x = evt[id]; // x, y, z...
            fptype offset = RO_CACHE(pc.parameters[pc.parameterIdx + 1 + j]); // x0, y0, z0...
            x -= offset;
            int currPower = currIndex % maxDegree;
            currIndex /= maxDegree;
            currTerm *= pow(x, currPower);
            sumOfIndices += currPower;
            // if ((gpuDebug & 1) && (THREADIDX == 50) && (BLOCKIDX == 3))
            // if ((BLOCKIDX == internalDebug1) && (THREADIDX == internalDebug2))
            // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
            // printf("  [%f -> %f^%i = %f] (%i %i) \n", evt[indices[2 + indices[0] + j]], x, currPower, pow(x,
            // currPower), sumOfIndices, indices[2 + indices[0] + j]);
        }

        // if ((gpuDebug & 1) && (THREADIDX == 50) && (BLOCKIDX == 3))
        // if ((BLOCKIDX == internalDebug1) && (THREADIDX == internalDebug2))
        // printf(") End box %i\n", i);
        // All threads should hit this at the same time and with the same result. No branching.
        if(sumOfIndices >= maxDegree)
            continue;

        fptype coefficient = RO_CACHE(pc.parameters[pc.parameterIdx + 1 + coeffNumber]); // Coefficient from MINUIT
        coeffNumber++;
        //if ((gpuDebug & 1) && (THREADIDX == 50) && (BLOCKIDX == 3))
        //if ((BLOCKIDX == internalDebug1) && (THREADIDX == internalDebug2))
        //if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
        //printf("Box %i contributes %f * %f = %f -> %f\n", i, currTerm, p[indices[coeffNumber - 1]], coefficient*currTerm, (ret + coefficient*currTerm));
        currTerm *= coefficient;
        ret += currTerm;
    }

    pc.incrementIndex (1, num_parameters, num_constants, num_observables, 1);

    return ret;
}

__device__ device_function_ptr ptr_to_Polynomial       = device_Polynomial;
__device__ device_function_ptr ptr_to_OffsetPolynomial = device_OffsetPolynomial;
__device__ device_function_ptr ptr_to_MultiPolynomial  = device_MultiPolynomial;

// Constructor for single-variate polynomial, with optional zero point.
__host__ PolynomialPdf::PolynomialPdf(
    std::string n, Variable *_x, std::vector<Variable *> weights, Variable *x0, unsigned int lowestDegree)
    : GooPdf(_x, n)
    , center(x0) {
    std::vector<unsigned int> pindices;
    pindices.push_back(lowestDegree);

    constantsList.push_back (lowestDegree);

    for(std::vector<Variable*>::iterator v = weights.begin(); v != weights.end(); ++v) {
        pindices.push_back(registerParameter(*v));
    }

    if(x0) {
        polyType = 1;
        pindices.push_back(registerParameter(x0));
        GET_FUNCTION_ADDR(ptr_to_OffsetPolynomial);
    } else {
        polyType = 0;
        GET_FUNCTION_ADDR(ptr_to_Polynomial);
    }

    initialize(pindices);
}

// Constructor for multivariate polynomial.
__host__ PolynomialPdf::PolynomialPdf(std::string n,
                                      std::vector<Variable *> obses,
                                      std::vector<Variable *> coeffs,
                                      std::vector<Variable *> offsets,
                                      unsigned int maxDegree)
    : GooPdf(nullptr, n) {
    unsigned int numParameters = 1;

    constantsList.push_back(maxDegree);

    // For 1 observable, equal to n = maxDegree + 1.
    // For two, n*(n+1)/2, ie triangular number. This generalises:
    // 3: Pyramidal number n*(n+1)*(n+2)/(3*2)
    // 4: Hyperpyramidal number n*(n+1)*(n+2)*(n+3)/(4*3*2)
    // ...
    for(unsigned int i = 0; i < obses.size(); ++i) {
        registerObservable(obses[i]);
        numParameters *= (maxDegree + 1 + i);

        //we are 'padding' the list.
        constantsList.push_back(maxDegree + 1 + i);
    }

    for(int i = observablesList.size(); i > 1; --i)
        numParameters /= i;

    while(numParameters > coeffs.size()) {
        char varName[100];
        sprintf(varName, "%s_extra_coeff_%i", getName().c_str(), static_cast<int>(coeffs.size()));

        Variable *newTerm = new Variable(varName, 0);
        coeffs.push_back(newTerm);

        std::cout << "Warning: " << getName() << " created dummy variable " << varName
                  << " (fixed at zero) to account for all terms.\n";
    }

    while(offsets.size() < obses.size()) {
        char varName[100];
        sprintf(varName, "%s_extra_offset_%i", getName().c_str(), static_cast<int>(offsets.size()));
        Variable *newOffset = new Variable(varName, 0);
        offsets.push_back(newOffset);
    }

    std::vector<unsigned int> pindices;
    pindices.push_back(maxDegree);

    for(auto &offset : offsets) {
        pindices.push_back(registerParameter(offset));
    }

    for(auto &coeff : coeffs) {
        pindices.push_back(registerParameter(coeff));
    }

    polyType = 2;
    GET_FUNCTION_ADDR(ptr_to_MultiPolynomial);
    initialize(pindices);
}

__host__ void PolynomialPdf::recursiveSetIndices () {
    if (polyType == 0)
    {
        GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName (), "ptr_to_Polynomial");
        GET_FUNCTION_ADDR(ptr_to_Polynomial);
    }
    else if (polyType == 1)
    {
        GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName (), "ptr_to_OffsetPolynomia");
        GET_FUNCTION_ADDR(ptr_to_OffsetPolynomial);
    }
    else if (polyType == 2)
    {
        GOOFIT_TRACE("host_function_table[{}] = {}({})", num_device_functions, getName (), "ptr_to_MultiPolynomial");
        GET_FUNCTION_ADDR(ptr_to_MultiPolynomial);
    }

    GOOFIT_TRACE("host_function_table[{}] = {}", num_device_functions, getName ());
    host_function_table[num_device_functions] = host_fcn_ptr;
    functionIdx = num_device_functions ++;

    populateArrays ();
}

__host__ fptype PolynomialPdf::integrate(fptype lo, fptype hi) const {
    // This is *still* wrong. (13 Feb 2013.)
    fptype lowestDegree = host_constants[constantsIdx + 1];

    if(center) {
        hi -= host_observables[observablesIdx + 1];
        lo -= host_observables[observablesIdx + 2];
    }

    fptype ret = 0;

    for(int i = 2; i < host_parameters[parametersIdx] + (center ? 0 : 1); ++i) {
        fptype powerPlusOne = lowestDegree + i - 2;
        fptype curr = pow(hi, powerPlusOne);
        curr       -= pow(lo, powerPlusOne);
        curr       /= powerPlusOne;
        ret        += host_parameters[parametersIdx + i] * curr;
    }

    return ret;
}

__host__ fptype PolynomialPdf::getCoefficient(int coef) const {
    // NB! This function only works for single polynomials.
    if(1 != observablesList.size()) {
        std::cout <<
                  "Warning: getCoefficient method of PolynomialPdf not implemented for multi-dimensional polynomials. Returning zero, which is very likely wrong.\n";
        return 0;
    }

    // True function is, say, ax^2 + bx + c.
    // We express this as (a'x^2 + b'x + c')*N.
    // So to get the true coefficient, multiply the internal
    // one by the normalisation. (In non-PDF cases the normalisation
    // equals one, which gives the same result.)

    // Structure is nP lowestdegree c1 c2 c3 nO o1
    if(coef < host_constants[constantsIdx + 1])
        return 0; // Less than least power.

    if(coef > host_constants[constantsIdx + 1] + (host_parameters[parametersIdx] - 1))
        return 0; // Greater than max power.

    fptype norm = normalize();
    norm        = (1.0 / norm);

    fptype param = host_parameters[parametersIdx + 2 + coef - int(host_constants[constantsIdx + 1])];
    return norm*param;
}
} // namespace GooFit
