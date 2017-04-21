#pragma once

#include "goofit/PDFs/GooPdf.h"

class BWPdf : public GooPdf {

public:
    BWPdf(std::string n, Variable* _x, Variable* m, Variable* s);
private:

};
