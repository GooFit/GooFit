/*
04/18/2016 Christoph Hasse

This file includes a converter from 16 values (4 4-momenta) to (the common set of) 5 parameters.
See UserUtils for a handy way to convert for tests.

*/

#pragma once

#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/GFunctional.h>
#include <mcbooster/GTypes.h>
#include <mcbooster/Vector4R.h>

namespace GooFit {

struct Dim5 : public mcbooster::IFunctionArray {
    Dim5() { dim = 4; }

    __host__ __device__ auto cosHELANG(const mcbooster::Vector4R p,
                                       const mcbooster::Vector4R q,
                                       const mcbooster::Vector4R d) -> mcbooster::GReal_t {
        mcbooster::GReal_t pd  = p * d;
        mcbooster::GReal_t pq  = p * q;
        mcbooster::GReal_t qd  = q * d;
        mcbooster::GReal_t mp2 = p.mass2();
        mcbooster::GReal_t mq2 = q.mass2();
        mcbooster::GReal_t md2 = d.mass2();

        return (pd * mq2 - pq * qd) / sqrt((pq * pq - mq2 * mp2) * (qd * qd - mq2 * md2));
    }

    __host__ __device__ auto phi(const mcbooster::Vector4R &p4_p,
                                 const mcbooster::Vector4R &p4_d1,
                                 const mcbooster::Vector4R &p4_d2,
                                 const mcbooster::Vector4R &p4_h1,
                                 const mcbooster::Vector4R &p4_h2) -> mcbooster::GReal_t {
        mcbooster::Vector4R p4_d1p, p4_h1p, p4_h2p, p4_d2p;

        mcbooster::Vector4R d1_perp, d1_prime, h1_perp;
        mcbooster::Vector4R D;

        D                      = p4_d1 + p4_d2;
        mcbooster::Vector4R D2 = p4_h1 + p4_h2;

        d1_perp = p4_d1 - (D.dot(p4_d1) / D.dot(D)) * D;
        h1_perp = p4_h1 - (D2.dot(p4_h1) / D2.dot(D2)) * D2;

        // orthogonal to both D and d1_perp

        d1_prime = D.cross(d1_perp);

        d1_perp  = d1_perp / d1_perp.d3mag();
        d1_prime = d1_prime / d1_prime.d3mag();

        mcbooster::GReal_t x, y;

        x                      = d1_perp.dot(h1_perp);
        y                      = d1_prime.dot(h1_perp);
        mcbooster::GReal_t phi = atan2(y, x);
        // printf("x:%.5g, y%.5g phi %.5g\n", x, y, phi );

        if(phi < 0.0)
            phi += 2.0 * M_PI;

        mcbooster::Vector4R d1n  = p4_d1 / p4_d1.d3mag();
        mcbooster::Vector4R d2n  = p4_d2 / p4_d2.d3mag();
        mcbooster::Vector4R h1n  = p4_h1 / p4_h1.d3mag();
        mcbooster::Vector4R h2n  = p4_h2 / p4_h2.d3mag();
        mcbooster::Vector4R h12n = (p4_h1 + p4_h2);
        h12n *= 1.0 / h12n.d3mag();

        mcbooster::Vector4R n1 = d1n.cross(d2n);
        mcbooster::Vector4R n2 = h1n.cross(h2n);
        n1 *= 1.0 / n1.d3mag();
        n2 *= 1.0 / n2.d3mag();
        mcbooster::Vector4R n3 = n1.cross(n2);

        mcbooster::GReal_t cp = (n1.dot(n2));
        mcbooster::GReal_t sp = (n3.dot(h12n));

        mcbooster::GReal_t phi2 = acos(cp);

        if(sp < 0)
            phi2 *= -1;

        // printf("cp %.5g, sp %.5g, phi2 %.5g\n",cp, sp, phi2 );

        return phi2;
    }

    __host__ __device__ void
    operator()(const mcbooster::GInt_t n, mcbooster::Vector4R **particles, mcbooster::GReal_t *variables) override {
        // int tid = blockIdx.x*blockDim.x + threadIdx.x;
        // bool PRINT_ME = tid == 0;
        const bool PRINT_ME = false;

        mcbooster::Vector4R ppip  = *particles[0];
        mcbooster::Vector4R ppim  = *particles[1];
        mcbooster::Vector4R pK    = *particles[2];
        mcbooster::Vector4R ppip2 = *particles[3];

        if (PRINT_ME)
        {
            printf("E, px, py, pz\n");
            printf("ppip = %f, %f, %f, %f\n", ppip.get(0), ppip.get(1), ppip.get(2), ppip.get(3));
            printf("ppim = %f, %f, %f, %f\n", ppim.get(0), ppim.get(1), ppim.get(2), ppim.get(3));
            printf("pK = %f, %f, %f, %f\n", pK.get(0), pK.get(1), pK.get(2), pK.get(3));
            printf("ppip2 = %f, %f, %f, %f\n", ppip2.get(0), ppip2.get(1), ppip2.get(2), ppip2.get(3));
        }

        mcbooster::Vector4R pM    = ppip + ppim + pK + ppip2;
        mcbooster::Vector4R ppipi = ppip + ppim;
        mcbooster::Vector4R pKpi  = pK + ppip2;

        if (PRINT_ME)
        {
            printf("pM = %f, %f, %f, %f\n", pM.get(0), pM.get(1), pM.get(2), pM.get(3));
            printf("ppipi = %f, %f, %f, %f\n", ppipi.get(0), ppipi.get(1), ppipi.get(2), ppipi.get(3));
            printf("pKpi = %f, %f, %f, %f\n", pKpi.get(0), pKpi.get(1), pKpi.get(2), pKpi.get(3));
        }

        variables[0] = ppipi.mass();
        variables[1] = pKpi.mass();
        variables[2] = cosHELANG(pM, ppipi, ppip);
        variables[3] = cosHELANG(pM, pKpi, pK);
        variables[4] = phi(pM, ppip, ppim, pK, ppip2);

        if (PRINT_ME)
        {
            printf("ppipi.mass() = %f\n", variables[0]);
            printf("pKpi.mass() = %f\n", variables[1]);
            printf("cosHELANG(pM, ppipi, ppip) = %f\n", variables[2]);
            printf("cosHELANG(pM, pKpi, pK) = %f\n", variables[3]);
            printf("phi = %f\n", variables[4]);
        }
    }
};

} // namespace GooFit
