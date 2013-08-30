#ifndef BINTRANSFORM_PDF_HH
#define BINTRANSFORM_PDF_HH

#include "EngineCore.hh" 

// Transforms ND coordinates into a single bin number. 
class BinTransformPdf : public EngineCore {
public:
  BinTransformPdf (std::string n, vector<Variable*> obses, vector<fptype> limits, vector<fptype> binSizes, vector<int> numBins); 

private:

};

#endif
