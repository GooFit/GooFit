#include "goofit/detail/CLI.hpp"

#ifdef GOOFIT_OMP
#include <omp.h>
#endif

class App : public CLI::App {

public:

    
    void parse(int argc, char** argv) {
        #ifdef GOOFIT_OMP
        MPI_Init(&argc, &argv);
        #endif

        CLI::App::parse(argc, argv);
    }

    ~App() {
        #ifdef GOOFIT_OMP
        MPI_Finalize();
        #endif

    }
}


