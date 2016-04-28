/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!
See *.cu file for more details
*/

#ifndef SPIN_FACTORS_HH
#define SPIN_FACTORS_HH

#include "DalitzPlotHelpers.hh"
#include "LineshapesPdf.hh"

typedef fptype (*spin_function_ptr) (fptype*, unsigned int*); 

enum class SF_4Body{
  DtoV1V2_V1toP1P2_V2toP3P4_S,
  DtoV1V2_V1toP1P2_V2toP3P4_P,
  DtoV1V2_V1toP1P2_V2toP3P4_D,
  DtoAP1_AtoVP2_VtoP3P4,
  DtoAP1_AtoVP2Dwave_VtoP3P4,
  DtoVS_VtoP1P2_StoP3P4,
  DtoV1P1_V1toV2P2_V2toP3P4,
  DtoAP1_AtoSP2_StoP3P4
};

class SpinFactor : public GooPdf {
  friend class DPPdf;

public:
  SpinFactor(std::string name, SF_4Body SF, unsigned int P0, unsigned int P1, unsigned int P2, unsigned int P3);
  void setConstantIndex (unsigned int idx) {host_indices[parameters + 1] = idx;}
  
private:
  unsigned int kind;
};


#endif
