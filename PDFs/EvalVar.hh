/*
04/18/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!
See *.cu file for more details
*/

#ifndef EVAL_VAR_HH
#define EVAL_VAR_HH

#include "DalitzPlotHelpers.hh" 
#include <mcbooster/GTypes.h>
#include <mcbooster/Vector4R.h>
#include <mcbooster/GContainers.h>
#include <mcbooster/GFunctional.h>


struct Dim5: public MCBooster::IFunctionArray
{
  Dim5()
  {
  dim = 4;
  }

  __host__   __device__ MCBooster::GReal_t cosHELANG(const MCBooster::Vector4R p, const MCBooster::Vector4R q,
    const MCBooster::Vector4R d)
  {
  MCBooster::GReal_t pd = p * d;
  MCBooster::GReal_t pq = p * q;
  MCBooster::GReal_t qd = q * d;
  MCBooster::GReal_t mp2 = p.mass2();
  MCBooster::GReal_t mq2 = q.mass2();
  MCBooster::GReal_t md2 = d.mass2();

  return (pd * mq2 - pq * qd)
      / sqrt((pq * pq - mq2 * mp2) * (qd * qd - mq2 * md2));

  }

  __host__   __device__ MCBooster::GReal_t phi(const MCBooster::Vector4R& p4_p,
    const MCBooster::Vector4R& p4_d1, const MCBooster::Vector4R& p4_d2, const MCBooster::Vector4R& p4_h1,
    const MCBooster::Vector4R& p4_h2)
  {

  MCBooster::Vector4R p4_d1p, p4_h1p, p4_h2p, p4_d2p;

  MCBooster::Vector4R d1_perp, d1_prime, h1_perp;
  MCBooster::Vector4R D;

  D = p4_d1 + p4_d2;
  MCBooster::Vector4R D2 = p4_h1 + p4_h2;

  d1_perp = p4_d1 - (D.dot(p4_d1) / D.dot(D)) * D;
  h1_perp = p4_h1 - (D2.dot(p4_h1) / D2.dot(D2)) * D2;

  // orthogonal to both D and d1_perp

  d1_prime = D.cross(d1_perp);

  d1_perp = d1_perp / d1_perp.d3mag();
  d1_prime = d1_prime / d1_prime.d3mag();

  MCBooster::GReal_t x, y;

  x = d1_perp.dot(h1_perp);
  y = d1_prime.dot(h1_perp);
  MCBooster::GReal_t phi = atan2(y, x);
  // printf("x:%.5g, y%.5g phi %.5g\n", x, y, phi );

  if (phi < 0.0)
    phi += 2.0 * CONST_PI;

  MCBooster::Vector4R d1n = p4_d1/p4_d1.d3mag();
  MCBooster::Vector4R d2n = p4_d2/p4_d2.d3mag();
  MCBooster::Vector4R h1n = p4_h1/p4_h1.d3mag();
  MCBooster::Vector4R h2n = p4_h2/p4_h2.d3mag();
  MCBooster::Vector4R h12n = (p4_h1+p4_h2);
  h12n*=1.0/h12n.d3mag();


  MCBooster::Vector4R n1 = d1n.cross(d2n);
  MCBooster::Vector4R n2 = h1n.cross(h2n);
  n1 *= 1.0/n1.d3mag();
  n2 *= 1.0/n2.d3mag();
  MCBooster::Vector4R n3 = n1.cross(n2);

  MCBooster::GReal_t cp = (n1.dot(n2));
  MCBooster::GReal_t sp = (n3.dot(h12n));

  MCBooster::GReal_t phi2 = acos(cp);    
  if (sp <0) phi2 *= -1;
  // printf("cp %.5g, sp %.5g, phi2 %.5g\n",cp, sp, phi2 );



  return phi2;

  }

  __host__ __device__
  void operator()(const MCBooster::GInt_t n, MCBooster::Vector4R** particles, MCBooster::GReal_t* variables)
  {
  MCBooster::Vector4R ppip = *particles[0];
  MCBooster::Vector4R ppim    = *particles[1];
  MCBooster::Vector4R pK   = *particles[2];
  MCBooster::Vector4R ppip2  = *particles[3];

  MCBooster::Vector4R pM = ppip + ppim + pK + ppip2;
  MCBooster::Vector4R ppipi = ppip + ppim;
  MCBooster::Vector4R pKpi = pK + ppip2;


  variables[0] = ppipi.mass();
  variables[1] = pKpi.mass();
  variables[2] = cosHELANG(pM, ppipi, ppip);
  variables[3] = cosHELANG(pM, pKpi, pK);
  variables[4] = phi(pM, ppip, ppim, pK, ppip2);
  }
};
 


#endif

