#include <goofit/PDFs/physics/lineshapes/GSpline.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

#include <Eigen/Core>
#include <Eigen/LU>

namespace GooFit {

__device__ __thrust_forceinline__ auto Q2(fptype Msq, fptype M1sq, fptype M2sq) -> fptype {
    return (Msq / 4. - (M1sq + M2sq) / 2. + (M1sq - M2sq) * (M1sq - M2sq) / (4 * Msq));
}

__device__ __thrust_forceinline__ auto BlattWeisskopf_Norm(const fptype z2, const fptype z02, unsigned int L)
    -> fptype {
    if(L == 0)
        return 1;
    else if(L == 1)
        return (1 + z02) / (1 + z2);
    else if(L == 2)
        return (z02 * z02 + 3 * z02 + 9) / (z2 * z2 + 3 * z2 + 9);
    else {
        abort(__FILE__, __LINE__, "Wrong value of L");
        return 0; // Can't reach
    }
}

__device__ auto getSpline(fptype x, bool continued, ParameterContainer &pc) -> fptype {
    const fptype s_min       = pc.getConstant(2);
    const fptype s_max       = pc.getConstant(3);
    const unsigned int nBins = pc.getConstant(4);
    const int numConstants   = pc.getNumConstants();

    constexpr int start_index = 5;

    // 11 is the first spine knot, 11+nBins is the first curvature

    if(x <= s_min) {
        fptype m_x_0 = pc.getConstant(start_index + 0);
        return continued ? m_x_0 : 0;
    } else if(x >= s_max) {
        fptype m_x_1 = pc.getConstant(start_index + nBins - 1);
        return continued ? m_x_1 : 0;
    }

    fptype spacing = (s_max - s_min) / (nBins - 1.);
    fptype dx      = fmod((x - s_min), spacing);

    auto bin = static_cast<unsigned int>((x - s_min) / spacing);

    fptype m_x_0  = pc.getConstant(start_index + bin);
    fptype m_x_1  = pc.getConstant(start_index + bin + 1);
    fptype m_xf_0 = pc.getConstant(start_index + bin + nBins);
    fptype m_xf_1 = pc.getConstant(start_index + bin + nBins + 1);

    // TODO: try to calculate this.
    pc.incrementIndex(1, 0, numConstants, 0, 1);

    return m_x_0 + dx * ((m_x_1 - m_x_0) / spacing - (m_xf_1 + 2 * m_xf_0) * spacing / 6) + dx * dx * m_xf_0
           + dx * dx * dx * (m_xf_1 - m_xf_0) / (6 * spacing);
}

__device__ auto kFactor(fptype mass, fptype width) -> fptype {
    fptype gamma = mass * sqrt(POW2(mass) + POW2(width));
    fptype k     = 2 * sqrt(2.) * mass * width * gamma / (M_PI * sqrt(POW2(mass) + gamma));
    return sqrt(k);
}

__device__ auto Spline_TDP(fptype Mpair, fptype m1, fptype m2, ParameterContainer &pc) -> fpcomplex {
    const fptype mass      = pc.getParameter(0);
    const fptype width     = pc.getParameter(1);
    const fptype radius    = pc.getConstant(1);
    const int numConstants = pc.getNumConstants();

    fptype s  = POW2(Mpair);
    fptype s1 = POW2(m1);
    fptype s2 = POW2(m2);

    // This is GSpline.EFF in AmpGen

    fptype q2 = fabs(Q2(s, s1, s2));

    // Non-EFF
    // fptype BF             = sqrt( BlattWeisskopf_Norm(q2 * POW2(radius), 0, L));
    fptype BF = exp(-q2 * POW2(radius) / 2);

    fptype width_shape = width * getSpline(s, true, pc);
    fptype width_norm  = width * getSpline(POW2(mass), false, pc);

    fptype norm          = kFactor(mass, width) * BF;
    fptype running_width = width_norm == 0 ? 0 : width * width_shape / width_norm;
    fpcomplex iBW        = fpcomplex(POW2(mass) - s, -mass * running_width);

    pc.incrementIndex(1, 2, numConstants, 0, 1);

    return norm / iBW;
}

__device__ resonance_function_ptr ptr_to_Spline = Spline_TDP;

auto make_spline_curvatures(std::vector<Variable> vars, Lineshapes::spline_t SplineInfo) -> std::vector<fptype> {
    size_t size = std::get<2>(SplineInfo) - 2;
    Eigen::Matrix<fptype, Eigen::Dynamic, Eigen::Dynamic> m(size, size);
    for(size_t i = 0; i < size; i++) {
        m(i, i) = 4;
        if(i != size - 1) {
            m(i, i + 1) = 1;
            m(i + 1, 1) = 1;
        }
    }
    m = m.inverse();

    Eigen::Matrix<fptype, Eigen::Dynamic, 1> L(size);
    for(unsigned int i = 0; i < size; ++i)
        L[i] = vars[i + 2].getValue() - 2 * vars[i + 1].getValue() + vars[i].getValue();

    auto mtv       = m * L;
    fptype spacing = (std::get<0>(SplineInfo) - std::get<1>(SplineInfo)) / std::get<2>(SplineInfo);

    std::vector<fptype> ret(vars.size(), 0);
    for(unsigned int i = 0; i < size; ++i) {
        ret.at(i + 1) = 6 * mtv(i) / POW2(spacing);
    }
    return ret;
}

Lineshapes::GSpline::GSpline(std::string name,
                             Variable mass,
                             Variable width,
                             unsigned int L,
                             unsigned int Mpair,
                             FF FormFac,
                             fptype radius,
                             std::vector<Variable> AdditionalVars,
                             spline_t SplineInfo)
    : Lineshape("GSpline", name, L, Mpair, FormFac, radius) {
    registerParameter(mass);
    registerParameter(width);

    registerConstant(radius);

    if(std::get<2>(SplineInfo) != AdditionalVars.size())
        throw GeneralError("bins {} != vars {}", std::get<2>(SplineInfo), AdditionalVars.size());
    registerConstant(std::get<0>(SplineInfo));
    registerConstant(std::get<1>(SplineInfo));
    registerConstant(std::get<2>(SplineInfo));

    for(auto &par : AdditionalVars) {
        registerParameter(par);
    }

    std::vector<fptype> SplineCTerms = make_spline_curvatures(AdditionalVars, SplineInfo);

    for(auto &par : SplineCTerms) {
        // Calc curve
        registerConstant(par);
    }

    registerFunction("ptr_to_Spline", ptr_to_Spline);

    initialize();
}

bool Lineshapes::GSpline::isEqualByValue(const Lineshape &other) const {
    return this->Lineshape::isEqualByValue(other);
}

} // namespace GooFit
