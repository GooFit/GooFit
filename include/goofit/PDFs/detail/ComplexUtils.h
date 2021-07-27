#pragma once

#include <goofit/Error.h>
#include <goofit/detail/Complex.h>
#include <vector>

namespace GooFit {

/// Flatten a complex array into a standard one (1r, 1i, 2r, 2i, ...)
template <typename T>
auto flatten(const std::vector<thrust::complex<T>> &input) -> std::vector<T> {
    std::vector<T> output;
    for(auto val : input) {
        output.push_back(val.real());
        output.push_back(val.imag());
    }
    return output;
}

auto complex_derivative(const std::vector<fptype> &x, const std::vector<fpcomplex> &y) -> std::vector<fpcomplex> {
    if(x.size() != y.size()) // Must be a valid pointer
        throw GeneralError("x and y must have the same diminsions!");

    int i, k;
    unsigned int n = x.size();
    std::vector<fpcomplex> u(n);
    std::vector<fpcomplex> y2(n);

    fptype sig, p, qn, un;
    fpcomplex yp1 = 2. * (y[1] - y[0]) / (x[1] - x[0]);
    fpcomplex ypn = 2. * (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);

    /* The lower boundary condition is set either to be "natural" or else to have specified first derivative*/
    if(yp1.real() > 0.99e30) {
        y2[0].real(0.);
        u[0].real(0.);
    } else {
        y2[0].real(-0.5);
        u[0].real(3.0 / (x[1] - x[0]) * ((y[1].real() - y[0].real()) / (x[1] - x[0]) - yp1.real()));
    }
    if(yp1.imag() > 0.99e30) {
        y2[0].imag(0.);
        u[0].imag(0.);
    } else {
        y2[0].imag(-0.5);
        u[0].imag(3.0 / (x[1] - x[0]) * ((y[1].imag() - y[0].imag()) / (x[1] - x[0]) - yp1.imag()));
    }

    /* This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the
     * decomposed factors*/

    for(i = 1; i < n - 1; i++) {
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        p   = sig * y2[i - 1].real() + 2.0;
        y2[i].real((sig - 1.0) / p);
        u[i].real((y[i + 1].real() - y[i].real()) / (x[i + 1] - x[i])
                  - (y[i].real() - y[i - 1].real()) / (x[i] - x[i - 1]));
        u[i].real((6.0 * u[i].real() / (x[i + 1] - x[i - 1]) - sig * u[i - 1].real()) / p);
        p = sig * y2[i - 1].imag() + 2.0;
        y2[i].imag((sig - 1.0) / p);
        u[i].imag((y[i + 1].imag() - y[i].imag()) / (x[i + 1] - x[i])
                  - (y[i].imag() - y[i - 1].imag()) / (x[i] - x[i - 1]));
        u[i].imag((6.0 * u[i].imag() / (x[i + 1] - x[i - 1]) - sig * u[i - 1].imag()) / p);
    }

    /* The upper boundary condition is set either to be "natural" or else to have specified first derivative*/

    if(ypn.real() > 0.99e30) {
        qn = 0.;
        un = 0.;
    } else {
        qn = 0.5;
        un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn.real() - (y[n - 1].real() - y[n - 2].real()) / (x[n - 1] - x[n - 2]));
    }
    y2[n - 1].real((un - qn * u[n - 2].real()) / (qn * y2[n - 2].real() + 1.0));
    if(ypn.imag() > 0.99e30) {
        qn = 0.;
        un = 0.;
    } else {
        qn = 0.5;
        un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn.imag() - (y[n - 1].imag() - y[n - 2].imag()) / (x[n - 1] - x[n - 2]));
    }
    y2[n - 1].imag((un - qn * u[n - 2].imag()) / (qn * y2[n - 2].imag() + 1.0));

    /* This is the backsubstitution loop of the tridiagonal algorithm */

    for(k = n - 2; k >= 0; k--) {
        y2[k].real(y2[k].real() * y2[k + 1].real() + u[k].real());
        y2[k].imag(y2[k].imag() * y2[k + 1].imag() + u[k].imag());
    }

    return y2;
}

} // namespace GooFit
