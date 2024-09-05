#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>
#include <goofit/PDFs/physics/resonances/ResonanceUtils.h>

namespace GooFit {

namespace Resonances {

/// Relativistic Breit-Wigner
class Rescattering2 : public ResonancePdf {
  public:
    Rescattering2(std::string name,
        Variable ar,
        Variable ai,
        std::vector<Variable> coefs,
        unsigned int cyc,
        bool sym = false);
    ~Rescattering2() override = default;
};

} // namespace Resonances

} // namespace GooFit

/* 
coefs follow this order:
    LauParameter* B1_; 
		LauParameter* B2_; 
		LauParameter* B3_; 
		LauParameter* C1_; 
		LauParameter* C2_; 
		LauParameter* C3_; 
		LauParameter* C4_; 
		LauParameter* C5_; 
		LauParameter* D0_; 
		LauParameter* D1_; 
		LauParameter* D2_; 
		LauParameter* D3_; 
		LauParameter* F1_; 
		LauParameter* F2_; 
		LauParameter* F3_; 
		LauParameter* F4_; 

default values:
    const Double_t B1Val(23.6);
    const Double_t B2Val(29.4);
    const Double_t B3Val(0.6);
    const Double_t C1Val(34.39);
    const Double_t C2Val(4.4);
    const Double_t C3Val(-32.9);
    const Double_t C4Val(-16.);
    const Double_t C5Val(7.4);
    const Double_t D0Val(0.59);
    const Double_t D1Val(-0.38);
    const Double_t D2Val(0.12);
    const Double_t D3Val(-0.09);
    const Double_t F1Val(-0.043);
    const Double_t F2Val(-0.008);
    const Double_t F3Val(-0.28);
    const Double_t F4Val(0.026);
*/