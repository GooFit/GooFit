#include "mpi.h"
#include "TestThrust.hh"
#include "TRandom.h"
#include <complex>
#include "Faddeeva.cc"
#include <fstream>
#include <cassert>
#include <cstdlib>

int rank = 0;
int nCPU = 1;
int length = 0;

double cpuvoigtian(double x, double m, double w, double s) {
    // This calculation includes the normalisation - integral
    // over the reals is equal to one.

    // return constant for zero width and sigma
    if((0==s) && (0==w))
        return 1;

    assert(s > 0);
    assert(w > 0);

    double coef = -0.5/(s*s);
    double arg = x - m;

    // Breit-Wigner for zero sigma
    if(0==s)
        return (1/(arg*arg+0.25*w*w));

    // Gauss for zero width
    if(0==w)
        return exp(coef*arg*arg);

    // actual Voigtian for non-trivial width and sigma
    double c = 1./(sqrt(2)*s);
    double a = 0.5*c*w;
    double u = c*arg;
    std::complex<double> z(u, a) ;
    //printf("Calling Faddeeva %f %f %f %f %f %f %f\n", x, m, s, w, c, a, u);
    std::complex<double> v = Faddeeva_2(z);

    static const double rsqrtPi = 0.5641895835477563;
    return c*rsqrtPi*v.real();
}

double generateVoig() {
    static TRandom donram(rank);
    static const double mean = 0;
    static const double width = 0.2;
    static const double sigma = 0.1;
    static const double maxval = (4/(width*width));

    while(true) {
        double candidate = mean + donram.Uniform()*(100*width)-(50*width);
        double roll = donram.Uniform()*maxval;
        double val = cpuvoigtian(candidate, mean, width, sigma);
        assert(val <= maxval);

        if(val < roll)
            continue;

        return candidate;
    }
}

int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    rank = MPI::COMM_WORLD.Get_rank();

    length = atoi(argv[1]);
    std::ofstream writer;
    char buffer[200];
    sprintf(buffer, "/data/rolfa/bwdata_%i", rank);
    writer.open(buffer);

    for(int i = 0; i < length; ++i) {
        double cand = generateVoig();
        writer << cand << " ";
        //if (0 == i % 1000000) std::cout << "Generated " << i << " events.\n";
    }

    writer.close();

    MPI::Finalize();
    return 0;
}
