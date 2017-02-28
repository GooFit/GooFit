#include "goofit/detail/CLI.hpp"

#ifdef GOOFIT_OMP
#include <omp.h>
#endif

#ifdef GOOFIT_MPI
#include <mpi.h>
#endif

namespace GooFit {

class Application : public CLI::App {

public:
    //using CLI::App::App;

    virtual void parse(int argc, char** argv) {
#ifdef GOOFIT_MPI
        MPI_Init(&argc, &argv);

	int myId, numProcs;
	MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank (MPI_COMM_WORLD, &myId);

#ifndef GOOFIT_OMP
	int deviceCount;
	cudaGetDeviceCount (&deviceCount);

	int nodes = atoi (getenv("PBS_NUM_NODES"));
	if (nodes == 0)
	    nodes = 1;

	int procsPerNode = numProcs/nodes;
	int localRank = myId % procsPerNode;

	//Note, this will (probably) be overwritten by gpu-set-device calls...
	if (deviceCount == 1 && localRank > 1)
	{
	    //Multi-process to one GPU!
	    cudaSetDevice(0);
	}
	else if (procsPerNode > 1 && deviceCount > 1)
	{
            if (localRank <= deviceCount)
	    {
		//setting multiple processes to multiple GPU's
		cudaSetDevice (localRank);
	    }
	    else
            {
		//More multiple processes than GPU's, distribute sort of evenly!
		cudaSetDevice (localRank % deviceCount);
	    }
	}
	else
	{
	    //multiple GPU's, using one process
	    cudaSetDevice(0);
	}
#endif
#endif

        CLI::App::parse(argc, argv);
    }

    Application (std::string desc) : CLI::App (desc) {
    }

    virtual ~Application() {
#ifdef GOOFIT_MPI
        MPI_Finalize();
#endif

    }
};

}
