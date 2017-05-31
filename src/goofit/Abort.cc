#include "goofit/detail/Abort.h"
#include "goofit/Color.h"
#include "goofit/Error.h"
#include "goofit/PdfBase.h"
#include "goofit/detail/Backtrace.h"

namespace GooFit {

void abort(std::string file, int line, std::string reason, const PdfBase *pdf) {
    void *stackarray[20];

    std::cout << GooFit::reset << GooFit::red << "Abort called from " << file << " line " << line << " due to "
              << reason << std::endl;

    if(pdf) {
        std::vector<Variable *> pars = pdf->getParameters();
        std::cout << "Parameters of " << pdf->getName() << " : \n";

        for(Variable *v : pars) {
            if(0 > v->getIndex())
                continue;

            std::cout << "  " << v->getName() << " (" << v->getIndex() << ") :\t" << host_params[v->getIndex()]
                      << std::endl;
        }
    }

    std::cout << "Parameters (" << totalParams << ") :\n";

    for(int i = 0; i < totalParams; ++i) {
        std::cout << host_params[i] << " ";
    }

#if Backtrace_FOUND
    std::cout << GooFit::bold << std::endl;
    // get void* pointers for all entries on the stack
    size_t size = backtrace(stackarray, 20);
    // print out all the frames to stderr
    backtrace_symbols_fd(stackarray, size, 2);
#endif

    std::cout << GooFit::reset << std::flush;

    throw GooFit::GeneralError(reason);
}

} // namespace GooFit
