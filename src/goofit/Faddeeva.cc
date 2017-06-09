#include "goofit/Faddeeva.h"
#include "goofit/Error.h"

#include <cmath>

#include <thrust/complex.h>

namespace GooFit {

#define M_2PI 6.28318530717958

static fptype n1[12] = {0.25, 1.0, 2.25, 4.0, 6.25, 9.0, 12.25, 16.0, 20.25, 25.0, 30.25, 36.0};
static fptype e1[12] = {0.7788007830714049,
                        0.3678794411714423,
                        1.053992245618643e-1,
                        1.831563888873418e-2,
                        1.930454136227709e-3,
                        1.234098040866795e-4,
                        4.785117392129009e-6,
                        1.125351747192591e-7,
                        1.605228055185612e-9,
                        1.388794386496402e-11,
                        7.287724095819692e-14,
                        2.319522830243569e-16};

// table 2: coefficients for h = 0.53
static fptype n2[12]
    = {0.2809, 1.1236, 2.5281, 4.4944, 7.0225, 10.1124, 13.7641, 17.9776, 22.7529, 28.09, 33.9889, 40.4496};
static fptype e2[12] = {0.7551038420890235,
                        0.3251072991205958,
                        7.981051630007964e-2,
                        1.117138143353082e-2,
                        0.891593719995219e-3,
                        4.057331392320188e-5,
                        1.052755021528803e-6,
                        1.557498087816203e-8,
                        1.313835773243312e-10,
                        6.319285885175346e-13,
                        1.733038792213266e-15,
                        2.709954036083074e-18};

// tables for Pade approximation
static fptype C[7] = {65536.0, -2885792.0, 69973904.0, -791494704.0, 8962513560.0, -32794651890.0, 175685635125.0};
static fptype D[7] = {192192.0, 8648640.0, 183783600.0, 2329725600.0, 18332414100.0, 84329104860.0, 175685635125.0};

thrust::complex<fptype> Faddeeva_2(const thrust::complex<fptype> &z) {
    fptype *n, *e, t, u, r, s, d, f, g, h;
    thrust::complex<fptype> c, d2, v, w;
    int i;

    s = norm(z); // Actually the square of the norm. Don't ask me, I didn't name the function.

    if(s < 1e-7) {
        // use Pade approximation
        thrust::complex<fptype> zz = z * z;
        v                          = exp(zz);
        c                          = C[0];
        d2                         = D[0];

        for(i = 1; i <= 6; i++) {
            c  = c * zz + C[i];
            d2 = d2 * zz + D[i];
        }

        w = fptype(1.0) / v + thrust::complex<fptype>(0.0, M_2_SQRTPI) * c / d2 * z * v;
        return w;
    }

    // use trapezoid rule
    // select default table 1
    n = n1;
    e = e1;
    r = M_1_PI * 0.5;

    // if z is too close to a pole select table 2
    if(fabs(z.imag()) < 0.01 && fabs(z.real()) < 6.01) {
        h = fabs(z.real()) * 2;
        // Equivalent to modf(h, &g). Do this way because nvcc only knows about double version of modf.
        g = floor(h);
        h -= g;

        if(h < 0.02 || h > 0.98) {
            n = n2;
            e = e2;
            r = M_1_PI * 0.53;
        }
    }

    d = (z.imag() - z.real()) * (z.imag() + z.real());
    f = 4 * z.real() * z.real() * z.imag() * z.imag();

    g = h = 0.0;

    for(i = 0; i < 12; i++) {
        t = d + n[i];
        u = e[i] / (t * t + f);
        g += (s + n[i]) * u;
        h += (s - n[i]) * u;
    }

    u = 1 / s;

    c = r * thrust::complex<fptype>(z.imag() * (u + 2.0 * g), z.real() * (u + 2.0 * h));

    if(z.imag() < M_2PI) {
        s = 2.0 / r;
        t = s * z.real();
        u = s * z.imag();
        s = sin(t);
        h = cos(t);
        f = exp(-u) - h;
        g = 2.0 * exp(d - u) / (s * s + f * f);
        u = 2.0 * z.real() * z.imag();
        h = cos(u);
        t = sin(u);
        c += g * thrust::complex<fptype>((h * f - t * s), -(h * s + t * f));
    }

    return c;
}

fptype cpuvoigtian(fptype x, fptype m, fptype w, fptype s) {
    // This calculation includes the normalisation - integral
    // over the reals is equal to one.

    // return constant for zero width and sigma
    if((0 == s) && (0 == w))
        return 1;

    if(s <= 0 || w <= 0)
        throw GooFit::GeneralError("s {} and w {} must be larger than 0", s, w);

    fptype coef = -0.5 / (s * s);
    fptype arg  = x - m;

    // Breit-Wigner for zero sigma
    if(0 == s)
        return (1 / (arg * arg + 0.25 * w * w));

    // Gauss for zero width
    if(0 == w)
        return exp(coef * arg * arg);

    // actual Voigtian for non-trivial width and sigma
    fptype c = 1. / (sqrt(2) * s);
    fptype a = 0.5 * c * w;
    fptype u = c * arg;
    thrust::complex<fptype> z(u, a);
    // printf("Calling Faddeeva %f %f %f %f %f %f %f\n", x, m, s, w, c, a, u);
    thrust::complex<fptype> v = Faddeeva_2(z);

    static const fptype rsqrtPi = 0.5641895835477563;
    return c * rsqrtPi * v.real();
}
} // namespace GooFit
