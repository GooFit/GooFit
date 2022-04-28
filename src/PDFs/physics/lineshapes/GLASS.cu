#include <goofit/PDFs/physics/lineshapes/GLASS.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

#include "Common.h"

namespace GooFit {

// generalized lass lineshape as implemented in MINT3 by Tim Evans. if F=R=1 and phiF=phiR=0 this is equal to normal
// lass as implemented in Mint3.
// The difference between this and lass mint is not quite clear to me. need to get back to this later.
__device__ auto glass_MINT3(fptype Mpair, fptype m1, fptype m2, ParameterContainer &pc) -> fpcomplex {
    unsigned int orbital = pc.getConstant(1);
    fptype meson_radius  = pc.getConstant(2);

    fptype resmass  = pc.getParameter(0);
    fptype reswidth = pc.getParameter(1);
    fptype a        = pc.getParameter(2);
    fptype r        = pc.getParameter(3);
    fptype phiF     = pc.getParameter(4);
    fptype phiR     = pc.getParameter(5);
    fptype F        = pc.getParameter(6);

    fptype rMass2 = Mpair * Mpair;

    fptype R = 1.0;

    fptype mpsq  = (m1 + m2) * (m1 + m2);
    fptype mmsq  = (m1 - m2) * (m1 - m2);
    fptype num   = (rMass2 - mpsq) * (rMass2 - mmsq);
    fptype num2  = (resmass * resmass - mpsq) * (resmass * resmass - mmsq);
    fptype pABSq = num / (4 * rMass2);
    fptype prSq  = fabs(num2 / (4 * resmass * resmass));

    fptype pratio = sqrt(pABSq / prSq);

    fptype pratio_to_2Jplus1 = 1;

    for(int i = 0; i < 2 * orbital + 1; i++) {
        pratio_to_2Jplus1 *= pratio;
    }

    fptype mratio = resmass / Mpair;
    fptype r2     = meson_radius * meson_radius;
    fptype thisFR = BL_PRIME(pABSq * r2, prSq * r2, orbital);
    fptype GofM   = reswidth * pratio_to_2Jplus1 * mratio * thisFR;

    fptype y          = 2.0 * a * sqrt(pABSq);
    fptype x          = 2.0 + a * r * pABSq;
    fptype scattphase = phiF + atan(y / x);
    fptype resphase   = phiR + atan(resmass * GofM / (resmass * resmass - rMass2));
    fptype rho        = 1.0 / sqrt(pABSq / rMass2);
    fpcomplex returnVal
        = (F * sin(scattphase) * fpcomplex(cos(scattphase), sin(scattphase))
           + R * sin(resphase) * fpcomplex(cos(resphase + 2 * scattphase), sin(resphase + 2 * scattphase)))
          * rho;

    pc.incrementIndex(1, 7, 3, 0, 1);

    return returnVal;
}

__device__ resonance_function_ptr ptr_to_glass3 = glass_MINT3;

Lineshapes::GLASS::GLASS(std::string name,
                         Variable mass,
                         Variable width,
                         unsigned int L,
                         unsigned int Mpair,
                         FF FormFac,
                         fptype radius,
                         std::vector<Variable> AdditionalVars)
    : Lineshape("GLASS", name, L, Mpair, FormFac, radius) {
    if(5 != AdditionalVars.size()) {
        throw GeneralError("It seems you forgot to provide the vector with the five necessary variables for GLASS, a, "
                           "r, phiF, phiR and F (in that order)");
    }

    registerParameter(mass);
    registerParameter(width);
    for(int i = 0; i < 5; i++)
        registerParameter(AdditionalVars[i]);

    registerConstant(L);
    registerConstant(radius);

    registerFunction("ptr_to_glass3", ptr_to_glass3);

    initialize();
}

std::ostream &Lineshapes::GLASS::print(std::ostream &out) const {
    std::string paramNames     = "";
    std::string lassParamNames = "";
    for(int p = 0; p < this->getParameters().size(); p++) {
        if(p > 1) {
            lassParamNames += this->getParameters()[p].getName();
            if(p < this->getParameters().size() - 1) {
                lassParamNames += ", ";
            }
        } else {
            paramNames += this->getParameters()[p].getName() + ", ";
        }
    }

    DP4Pair mpairEnum = static_cast<DP4Pair>(this->_Mpair);

    out << this->getPdfName() << ", " << this->getName() << ", " << paramNames << this->_L << ", " << mpairEnum << ", "
        << this->_FormFac << ", " << this->_radius << ", " << lassParamNames;
    return out;
}

bool Lineshapes::GLASS::isEqualByValue(const Lineshape &other) const { return this->Lineshape::isEqualByValue(other); }

} // namespace GooFit
