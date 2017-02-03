#include "goofit/detail/CLI.hpp"

#ifdef GOOFIT_OMP
#include <omp.h>
#endif

namespace GooFit {

class Application : public CLI::App {

public:
    using CLI::App::App;

    void parse(int argc, char** argv) {
        #ifdef GOOFIT_OMP
        MPI_Init(&argc, &argv);
        #endif

        CLI::App::parse(argc, argv);
    }

    ~Application() {
        #ifdef GOOFIT_OMP
        MPI_Finalize();
        #endif

    }
};

}
