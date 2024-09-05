#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>
#include <goofit/PDFs/physics/resonances/ResonanceUtils.h>

namespace GooFit {

namespace Resonances {

/// POLE
		class Pole : public ResonancePdf {
			public:
				Pole(std::string name,
						Variable ar,
						Variable ai,
						Variable real,
						Variable img,
						unsigned int cyc,
						bool symmDP = false);
				~Pole() override = default;

    } ;

}// namespace Resonances

} // namespace GooFit
