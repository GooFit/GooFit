#pragma once

#include "goofit/PDFs/GooPdf.h"

class ExpGausPdf : public GooPdf {
public:
    ExpGausPdf(std::string n, Variable* _x, Variable* m, Variable* s, Variable* t);

private:

};
