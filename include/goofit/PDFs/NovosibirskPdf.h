#pragma once

#include "goofit/PDFs/GooPdf.h"

class NovosibirskPdf : public GooPdf {
public:
    NovosibirskPdf(std::string n, Variable* _x, Variable* m, Variable* s, Variable* t);

private:

};
