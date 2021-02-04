#pragma once

#include <goofit/PDFs/GooPdf.h>

namespace GooFit {

class PoissonPdf : public GooPdf {
	public:
		PoissonPdf(std::string n, Observable _x, Variable lambda);
};

}
