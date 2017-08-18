#include "goofit/PDFs/basic/PolynomialPdf.h"
#include "goofit/Variable.h"

namespace GooFit {

__device__ fptype device_Polynomial(fptype *evt, fptype *p, unsigned int *indices) {
    // Structure is nP lowestdegree c1 c2 c3 nO o1

    int numParams    = RO_CACHE(indices[0]) + 1;
    int lowestDegree = RO_CACHE(indices[1]);

    fptype x   = evt[RO_CACHE(indices[2 + RO_CACHE(indices[0])])];
    fptype ret = 0;

    for(int i = 2; i < numParams; ++i) {
        ret += RO_CACHE(p[RO_CACHE(indices[i])]) * pow(x, lowestDegree + i - 2);
    }

    return ret;
}

__device__ fptype device_OffsetPolynomial(fptype *evt, fptype *p, unsigned int *indices) {
    int numParams    = RO_CACHE(indices[0]);
    int lowestDegree = RO_CACHE(indices[1]);

    fptype x = evt[RO_CACHE(indices[2 + numParams])];
    x -= RO_CACHE(p[RO_CACHE(indices[numParams])]);
    fptype ret = 0;

    for(int i = 2; i < numParams; ++i) {
        ret += RO_CACHE(p[RO_CACHE(indices[i])]) * pow(x, lowestDegree + i - 2);
    }

    return ret;
}

__device__ fptype device_MultiPolynomial(fptype *evt, fptype *p, unsigned int *indices) {
    // Structure is nP, maxDegree, offset1, offset2, ..., coeff1, coeff2, ..., nO, o1, o2, ...

    int numObservables = RO_CACHE(indices[RO_CACHE(indices[0]) + 1]);
    int maxDegree      = RO_CACHE(indices[1]) + 1;
    // Only appears in construction (maxDegree + 1) or (x > maxDegree), so
    // may as well add the one and use >= instead.

    // Technique is to iterate over the full n-dimensional box, skipping matrix elements
    // whose sum of indices is greater than maxDegree. Notice that this is increasingly
    // inefficient as n grows, since a larger proportion of boxes will be skipped.
    int numBoxes = 1;

    for(int i = 0; i < numObservables; ++i)
        numBoxes *= maxDegree;

    int coeffNumber = 2 + numObservables; // Index of first coefficient is 2 + nO, not 1 + nO, due to maxDegree. (nO
                                          // comes from offsets.)
    fptype ret = RO_CACHE(p[RO_CACHE(indices[coeffNumber++])]); // Coefficient of constant term.

    for(int i = 1; i < numBoxes;
        ++i) { // Notice skip of inmost 'box' in the pyramid, corresponding to all powers zero, already accounted for.
        fptype currTerm  = 1;
        int currIndex    = i;
        int sumOfIndices = 0;

        // if ((gpuDebug & 1) && (THREADIDX == 50) && (BLOCKIDX == 3))
        // if ((BLOCKIDX == internalDebug1) && (THREADIDX == internalDebug2))
        // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
        // printf("[%i, %i] Start box %i %f %f:\n", BLOCKIDX, THREADIDX, i, ret, evt[8]);
        for(int j = 0; j < numObservables; ++j) {
            fptype x      = RO_CACHE(evt[RO_CACHE(indices[2 + RO_CACHE(indices[0]) + j])]); // x, y, z...
            fptype offset = RO_CACHE(p[RO_CACHE(indices[2 + j])]);                          // x0, y0, z0...
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

        fptype coefficient = RO_CACHE(p[RO_CACHE(indices[coeffNumber++])]); // Coefficient from MINUIT
        // if ((gpuDebug & 1) && (THREADIDX == 50) && (BLOCKIDX == 3))
        // if ((BLOCKIDX == internalDebug1) && (THREADIDX == internalDebug2))
        // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
        // printf("Box %i contributes %f * %f = %f -> %f\n", i, currTerm, p[indices[coeffNumber - 1]],
        // coefficient*currTerm, (ret + coefficient*currTerm));
        currTerm *= coefficient;
        ret += currTerm;
    }

    // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
    // printf("Final polynomial: %f\n", ret);

    // if (0 > ret) ret = 0;
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

    for(auto &weight : weights) {
        pindices.push_back(registerParameter(weight));
    }

    if(x0) {
        pindices.push_back(registerParameter(x0));
        GET_FUNCTION_ADDR(ptr_to_OffsetPolynomial);
    } else {
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

    // For 1 observable, equal to n = maxDegree + 1.
    // For two, n*(n+1)/2, ie triangular number. This generalises:
    // 3: Pyramidal number n*(n+1)*(n+2)/(3*2)
    // 4: Hyperpyramidal number n*(n+1)*(n+2)*(n+3)/(4*3*2)
    // ...
    for(unsigned int i = 0; i < obses.size(); ++i) {
        registerObservable(obses[i]);
        numParameters *= (maxDegree + 1 + i);
    }

    for(int i = observables.size(); i > 1; --i)
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

    GET_FUNCTION_ADDR(ptr_to_MultiPolynomial);
    initialize(pindices);
}

__host__ fptype PolynomialPdf::integrate(fptype lo, fptype hi) const {
    // This is *still* wrong. (13 Feb 2013.)
    unsigned int *indices = host_indices + parameters;
    fptype lowestDegree   = indices[1];

    if(center) {
        hi -= host_params[indices[indices[0]]];
        lo -= host_params[indices[indices[0]]];
    }

    fptype ret = 0;

    for(int i = 2; i < indices[0] + (center ? 0 : 1); ++i) {
        fptype powerPlusOne = lowestDegree + i - 2;
        fptype curr         = pow(hi, powerPlusOne);
        curr -= pow(lo, powerPlusOne);
        curr /= powerPlusOne;
        ret += host_params[indices[i]] * curr;
    }

    return ret;
}

__host__ fptype PolynomialPdf::getCoefficient(int coef) const {
    // NB! This function only works for single polynomials.
    if(1 != observables.size()) {
        std::cout << "Warning: getCoefficient method of PolynomialPdf not implemented for multi-dimensional "
                     "polynomials. Returning zero, which is very likely wrong.\n";
        return 0;
    }

    unsigned int *indices = host_indices + parameters;

    // True function is, say, ax^2 + bx + c.
    // We express this as (a'x^2 + b'x + c')*N.
    // So to get the true coefficient, multiply the internal
    // one by the normalisation. (In non-PDF cases the normalisation
    // equals one, which gives the same result.)

    // Structure is nP lowestdegree c1 c2 c3 nO o1
    if(coef < indices[1])
        return 0; // Less than least power.

    if(coef > indices[1] + (indices[0] - 1))
        return 0; // Greater than max power.

    fptype norm = normalize();
    norm        = (1.0 / norm);

    fptype param = host_params[indices[2 + coef - indices[1]]];
    return norm * param;
}
} // namespace GooFit
