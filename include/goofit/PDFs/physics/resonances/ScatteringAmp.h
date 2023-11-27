#pragma once
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

namespace Resonances {

class ScatteringAmp : public ResonancePdf {
			public:
				ScatteringAmp(std::string name,
						Variable ar,
						Variable ai,
						Variable akk,
						Variable apipi,
						std::vector<std::pair<fptype,fpcomplex>> &_phi00, // \delta_{KK}
						bool charge_pos,
						unsigned int cyc,
						bool symmDP = false);
				~ScatteringAmp() override = default;


};

}
}