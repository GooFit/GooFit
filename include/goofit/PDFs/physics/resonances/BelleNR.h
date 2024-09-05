#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>
#include <goofit/PDFs/physics/resonances/ResonanceUtils.h>

namespace GooFit {

namespace Resonances {

/// BELLE NR
		class BelleNR : public ResonancePdf {
			public:
				BelleNR(std::string name,
						Variable ar,
						Variable ai,
						Variable alpha,
						unsigned int cyc,
						bool symmDP = false);
				~BelleNR() override = default;

    } ;

}// namespace Resonances

} // namespace GooFit
