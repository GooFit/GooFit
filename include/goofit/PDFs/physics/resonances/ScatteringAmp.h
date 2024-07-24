#pragma once
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

namespace Resonances {

class ScatteringAmp : public ResonancePdf {
			public:
				ScatteringAmp(std::string name,
						Variable ar,
						Variable ai,
						std::complex<fptype> akk,
						std::complex<fptype> apipi,
						std::complex<fptype> f0scale,
						std::vector<std::pair<fptype,std::complex<fptype>>> _phi00, // \delta_{KK}
						unsigned int charge_pos,
						unsigned int cyc,
						bool symmDP = false);
				~ScatteringAmp() override = default;


};

}
}