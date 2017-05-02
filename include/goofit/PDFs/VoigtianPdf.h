#pragma once

#include "goofit/PDFs/GooPdf.h"

class VoigtianPdf : public GooPdf {
public:
    VoigtianPdf(std::string n, Variable* _x, Variable* m, Variable* s, Variable* w);

private:

};
