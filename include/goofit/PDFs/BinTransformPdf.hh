#ifndef BINTRANSFORM_PDF_HH
#define BINTRANSFORM_PDF_HH

#include "GooPdf.hh" 

// Transforms ND coordinates into a single bin number. 
class BinTransformPdf : public GooPdf {
public:
  BinTransformPdf (std::string n, vector<Variable*> obses, vector<fptype> limits, vector<fptype> binSizes, vector<int> numBins); 

private:

};

#endif
