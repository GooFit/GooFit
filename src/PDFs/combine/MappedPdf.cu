#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/combine/MappedPdf.h>

namespace GooFit {

__device__ fptype device_Mapped(fptype *evt, ParameterContainer &pc) {
    // Structure : nP mapFunctionIndex mapParamIndex functionIndex1 parameterIndex1 functionIndex2 parameterIndex2 ...

    // Find mapping between event variables and function to evaluate
    unsigned int numTargets = pc.getConstant(0);
    auto numConstants = pc.getNumConstants();

    auto pc_mapped = pc;

    // Mapping PDF happens directly after, so just increment.
    // pc.incrementIndex (1, 0, 1, 0, 1);
    pc.incrementIndex();

    // This is an index into the MappedPdf's list of functions
    // int targetFunction = (int) floor(0.5 +
    // (*(reinterpret_cast<device_function_ptr>(d_function_table[mapFunction])))(evt, p, paramIndices +
    // indices[2]));
    auto targetFunction = static_cast<int>(floor(0.5 + callFunction(evt, pc)));
    int indicestoskip = 0;
    for(int i = 0; i < targetFunction; i++ ) {
        indicestoskip += (int) pc_mapped.getConstant(1 + i);
    }
    //printf("indices to skip %d %d \n", targetFunction, indicestoskip );
   
    auto cur_funcIdx = pc.funcIdx;
    auto target_funcIdx = cur_funcIdx + indicestoskip;
    while(pc.funcIdx < target_funcIdx) pc.incrementIndex();
    //printf("index of fucntion %d \n", pc.funcIdx); 


    // targetFunction *= 2; // Because there are two pieces of information about each function
    // targetFunction += 3; // Because first function information begins at index 3

    // Unsure these carry over...
    // targetFunction *= 2;
    // targetFunction += 3;

    // numTargets; increment past our set of mapping functions and our target functions after it is handled.
    unsigned int funcIdx = pc.funcIdx;

    /*
    if(targetFunction == 0) pc.funcIdx = 2;
    if(targetFunction == 1) pc.funcIdx = 20;
    // brute force to make it work
    if(targetFunction == 1) {
        pc.parameterIdx += 154 - 2;
        pc.constantIdx += 63 - 6;
        pc.observableIdx += 38 - 8;
        pc.normalIdx += 40 - 4;
    }
    */


    // fptype ret = (*(reinterpret_cast<device_function_ptr>(d_function_table[indices[targetFunction]])))(evt, p,
    // paramIndices + indices[targetFunction + 1]);
    fptype norm = pc.getNormalization(0);
    fptype ret  = callFunction(evt, pc);
    //std::cout << "mapped return: " << norm << " " << ret << std::endl;
    ret *= norm;

    // increment our functions here...
    /*
    if (pc.funcIdx==20) pc.funcIdx=38;
    if (pc.funcIdx==19) pc.funcIdx=38;
    //DOES THIS HELP?
    if(targetFunction == 0) {
        pc.parameterIdx = 306;
        pc.constantIdx = 120;
        pc.observableIdx = 68;
        pc.normalIdx = 76;
    }
    */
    int finalIndex = cur_funcIdx;
    //now need to increase index until end
    for(int i = 1; i < numConstants; i++ ) {
        finalIndex += (int) pc_mapped.getConstant(i);
    }

    while(pc.funcIdx < finalIndex) pc.incrementIndex();
    //printf("final index %d\n", pc.funcIdx);



    // if (gpuDebug & 1)
    // if ((gpuDebug & 1) && (0 == BLOCKIDX) && (0 == THREADIDX))
    // printf("[%i, %i] Mapped: %i (%f %f %f %f) %f\n", BLOCKIDX, THREADIDX, targetFunction, evt[0], evt[1], evt[2],
    // evt[3], ret);
    return ret;
}

__device__ device_function_ptr ptr_to_Mapped = device_Mapped;

__host__ int countComponents(PdfBase *func) {
    auto subcomponents = func->getComponents();
    int n_components = 0;
    if(subcomponents.size() > 0) {
      
      for(auto subcomponent : subcomponents) {
         n_components++;
         //printf("%s\n", subcomponent->getName().c_str());
        n_components += countComponents(subcomponent);
      }
    }
    
    return n_components;
}

__host__ MappedPdf::MappedPdf(std::string n, GooPdf *m, std::vector<GooPdf *> &t)
    : CombinePdf("MappedPdf", n) {
    components.push_back(m);

    std::vector<int> nComponents;

    for(GooPdf *f : t) {
        components.push_back(f);
        // pindices.push_back(f->getFunctionIndex());
        // pindices.push_back(f->getParameterIndex());
        // functionIndicesUsed.insert(f->getFunctionIndex());
        auto functionIndex = f->getFunctionIndex();

        //printf("function index %d\n", functionIndex);
        //for(auto subcomp : f->getComponents()) printf("subcomp\n");
        //count number of subfunctions
        int n_components = countComponents(f);
        //also count total function
        n_components++;
        nComponents.push_back(n_components);
        //printf("nubmer comp %d\n", n_components);
    }

    // if(functionIndicesUsed.size() > 1) {
    //    std::cout << "Warning: More than one function type given to MappedPdf " << getName()
    //              << " constructor. This may slow execution by causing sequential evaluations.\n";
    //}

    // This makes sure we have the appropriate amount of obs in our structure
    observablesList = getObservables();

    // add a constant value for the number of 't' functions, skipping 'm'.
    registerConstant(components.size() - 1);
    for(auto nComponent : nComponents) registerConstant(nComponent);

    registerFunction("ptr_to_Mapped", ptr_to_Mapped);

    initialize();
}

__host__ fptype MappedPdf::normalize() {
    // std::cout << "Normalizing MappedPdf " << getName() << std::endl;
    fptype ret = 0;

    for(unsigned int i = 1; i < components.size(); ++i) { // No need to normalize mapping function.
        fptype curr = components[i]->normalize();
        ret += curr;
    }

    host_normalizations[normalIdx + 1] = 1.0;
    cachedNormalization                = 1.0;

    return ret;
}
} // namespace GooFit
