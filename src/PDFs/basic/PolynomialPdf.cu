#include <goofit/Log.h>
#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/Variable.h>

namespace GooFit {

__device__ auto device_Polynomial(fptype *evt, ParameterContainer &pc) -> fptype {
    int id = pc.getObservable(0);
    // Structure is nP lowestdegree c1 c2 c3 nO o1

    int numParams    = pc.getNumParameters();
    int lowestDegree = pc.getConstant(0);

    fptype x   = RO_CACHE(evt[id]);
    fptype ret = 0;

    // unsure why this starts at i=2...
    for(int i = 0; i < numParams; ++i) {
        fptype param = pc.getParameter(i);
        ret += param * pow(x, lowestDegree + i);
    }

    pc.incrementIndex(1, numParams, 1, 1, 1);
    return ret;
}

__device__ auto device_OffsetPolynomial(fptype *evt, ParameterContainer &pc) -> fptype {
    int id = pc.getObservable(0);

    int numParams    = pc.getNumParameters();
    int lowestDegree = pc.getConstant(0);

    fptype x = RO_CACHE(evt[id]);
    // TODO: Not sure where this is pointing...
    // x -= RO_CACHE(p[RO_CACHE(indices[numParams])]);
    fptype ret = 0;

    for(int i = 2; i < numParams; ++i) {
        ret += pc.getParameter(i) * pow(x, lowestDegree + i - 2);
    }

    pc.incrementIndex(1, numParams, 1, 1, 1);
    return ret;
}

__device__ auto device_MultiPolynomial(fptype *evt, ParameterContainer &pc) -> fptype {
    int num_constants  = pc.getNumConstants();
    int num_parameters = pc.getNumParameters();
    // Structure is nP, maxDegree, offset1, offset2, ..., coeff1, coeff2, ..., nO, o1, o2, ...

    int num_observables = pc.getNumObservables();
    int maxDegree       = pc.getConstant(0) + 1;
    // Only appears in construction (maxDegree + 1) or (x > maxDegree), so
    // may as well add the one and use >= instead.

    // Technique is to iterate over the full n-dimensional box, skipping matrix elements
    // whose sum of indices is greater than maxDegree. Notice that this is increasingly
    // inefficient as n grows, since a larger proportion of boxes will be skipped.
    int numBoxes = 1;

    for(int i = 0; i < num_observables; ++i)
        numBoxes *= maxDegree;

    int coeffNumber = num_observables; // Index of first coefficient is 2 + nO, not 1 + nO, due to maxDegree. (nO comes
                                       // from offsets.)
    fptype ret = pc.getParameter(coeffNumber); // Coefficient of constant term.
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
            // TODO:Need to debug these
            int id        = pc.getObservable(j);
            fptype x      = RO_CACHE(evt[id]);  // x, y, z...
            fptype offset = pc.getParameter(j); // x0, y0, z0...
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

        fptype coefficient = pc.getParameter(coeffNumber); // Coefficient from MINUIT
        coeffNumber++;
        // if ((gpuDebug & 1) && (THREADIDX == 50) && (BLOCKIDX == 3))
        // if ((BLOCKIDX == internalDebug1) && (THREADIDX == internalDebug2))
        // if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
        // printf("Box %i contributes %f * %f = %f -> %f\n", i, currTerm, p[indices[coeffNumber - 1]],
        // coefficient*currTerm, (ret + coefficient*currTerm));
        currTerm *= coefficient;
        ret += currTerm;
    }

    pc.incrementIndex(1, num_parameters, num_constants, num_observables, 1);

    return ret;
}

__device__ device_function_ptr ptr_to_Polynomial       = device_Polynomial;
__device__ device_function_ptr ptr_to_OffsetPolynomial = device_OffsetPolynomial;
__device__ device_function_ptr ptr_to_MultiPolynomial  = device_MultiPolynomial;

// Constructor for single-variate polynomial, with optional zero point.
__host__
PolynomialPdf::PolynomialPdf(std::string n, Observable _x, std::vector<Variable> weights, unsigned int lowestDegree)
    : GooPdf("PolynomialPdf", n, _x) {
    registerConstant(lowestDegree);

    for(Variable &v : weights) {
        registerParameter(v);
    }

    registerFunction("ptr_to_Polynomial", ptr_to_Polynomial);

    initialize();
}

__host__ PolynomialPdf::PolynomialPdf(
    std::string n, Observable _x, std::vector<Variable> weights, Variable x0, unsigned int lowestDegree)
    : GooPdf("PolynomialPdf", n, _x)
    , center(new Variable(x0)) {
    registerConstant(lowestDegree);

    for(Variable &v : weights) {
        registerParameter(v);
    }

    registerParameter(x0);

    registerFunction("ptr_to_OffsetPolynomial", ptr_to_OffsetPolynomial);

    initialize();
}

// Constructor for multivariate polynomial.
__host__ PolynomialPdf::PolynomialPdf(std::string n,
                                      std::vector<Observable> obses,
                                      std::vector<Variable> coeffs,
                                      std::vector<Variable> offsets,
                                      unsigned int maxDegree)
    : GooPdf("PolynomialPdf", n) {
    unsigned int numParameters = 1;

    registerConstant(maxDegree);
    // For 1 observable, equal to n = maxDegree + 1.
    // For two, n*(n+1)/2, ie triangular number. This generalises:
    // 3: Pyramidal number n*(n+1)*(n+2)/(3*2)
    // 4: Hyperpyramidal number n*(n+1)*(n+2)*(n+3)/(4*3*2)
    // ...
    for(unsigned int i = 0; i < obses.size(); ++i) {
        registerObservable(obses[i]);
        numParameters *= (maxDegree + 1 + i);
        // we are 'padding' the list.
        registerConstant(maxDegree + 1 + i);
    }

    for(int i = observablesList.size(); i > 1; --i)
        numParameters /= i;

    while(numParameters > coeffs.size()) {
        char varName[100];
        sprintf(varName, "%s_extra_coeff_%i", getName().c_str(), static_cast<int>(coeffs.size()));

        Variable newTerm(varName, 0);
        coeffs.push_back(newTerm);

        std::cout << "Warning: " << getName() << " created dummy variable " << varName
                  << " (fixed at zero) to account for all terms.\n";
    }

    while(offsets.size() < obses.size()) {
        char varName[100];
        sprintf(varName, "%s_extra_offset_%i", getName().c_str(), static_cast<int>(offsets.size()));
        Variable newOffset(varName, 0);
        offsets.push_back(newOffset);
    }

    for(auto &offset : offsets) {
        registerParameter(offset);
    }

    for(auto &coeff : coeffs) {
        registerParameter(coeff);
    }

    registerFunction("ptr_to_MultiPolynomial", ptr_to_MultiPolynomial);

    initialize();
}

__host__ auto PolynomialPdf::integrate(fptype lo, fptype hi) const -> fptype {
    // This is *still* wrong. (13 Feb 2013.)
    fptype lowestDegree = host_constants[constantsIdx + 1];

    if(center) {
        hi -= host_observables[observablesIdx + 1];
        lo -= host_observables[observablesIdx + 2];
    }

    fptype ret = 0;

    for(int i = 2; i < host_parameters[parametersIdx] + (center ? 0 : 1); ++i) {
        fptype powerPlusOne = lowestDegree + i - 2;
        fptype curr         = pow(hi, powerPlusOne);
        curr -= pow(lo, powerPlusOne);
        curr /= powerPlusOne;
        ret += host_parameters[parametersIdx + i] * curr;
    }

    return ret;
}

__host__ auto PolynomialPdf::getCoefficient(int coef) -> fptype {
    // NB! This function only works for single polynomials.
    if(1 != observablesList.size()) {
        std::cout << "Warning: getCoefficient method of PolynomialPdf not implemented for multi-dimensional "
                     "polynomials. Returning zero, which is very likely wrong.\n";
        return 0;
    }

    // True function is, say, ax^2 + bx + c.
    // We express this as (a'x^2 + b'x + c')*N.
    // So to get the true coefficient, multiply the internal
    // one by the normalization. (In non-PDF cases the normalization
    // equals one, which gives the same result.)

    // Structure is nP lowestdegree c1 c2 c3 nO o1
    if(coef < host_constants[constantsIdx + 1])
        return 0; // Less than least power.

    if(coef > host_constants[constantsIdx + 1] + (host_parameters[parametersIdx] - 1))
        return 0; // Greater than max power.

    fptype norm = normalize();
    norm        = (1.0 / norm);

    fptype param = host_parameters[parametersIdx + 2 + coef - int(host_constants[constantsIdx + 1])];
    return norm * param;
}
} // namespace GooFit
