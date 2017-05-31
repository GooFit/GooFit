/*
04/05/2016 Christoph Hasse
DISCLAIMER:

This code is not sufficently tested yet and still under heavy development!
See *.cu file for more details
*/

#pragma once

#include "goofit/PDFs/physics/DalitzPlotHelpers.h"

namespace GooFit {

typedef fptype (*spin_function_ptr)(fptype *, unsigned int *);

enum class SF_4Body {
    DtoPP1_PtoSP2_StoP3P4,
    DtoPP1_PtoVP2_VtoP3P4,
    DtoV1V2_V1toP1P2_V2toP3P4_S,
    DtoV1V2_V1toP1P2_V2toP3P4_P,
    DtoV1V2_V1toP1P2_V2toP3P4_D,
    DtoAP1_AtoVP2_VtoP3P4,
    DtoAP1_AtoVP2Dwave_VtoP3P4,
    DtoVS_VtoP1P2_StoP3P4,
    DtoV1P1_V1toV2P2_V2toP3P4,
    DtoAP1_AtoSP2_StoP3P4,
    DtoTP1_TtoVP2_VtoP3P4,
    FF_12_34_L1,
    FF_12_34_L2,
    FF_123_4_L1,
    FF_123_4_L2,
    ONE
};

class SpinFactor : public GooPdf {
    friend class DPPdf;

  private:
    SF_4Body _SF;
    unsigned int _P0;
    unsigned int _P1;
    unsigned int _P2;
    unsigned int _P3;

  public:
    SpinFactor(std::string name, SF_4Body SF, unsigned int P0, unsigned int P1, unsigned int P2, unsigned int P3);
    void setConstantIndex(unsigned int idx) { host_indices[parameters + 1] = idx; }
    bool operator==(const SpinFactor &S) const {
        return (S.getName() == getName() and S._SF == _SF and S._P0 == _P0 and S._P1 == _P1 and S._P2 == _P2
                and S._P3 == _P3);
    }
};
} // namespace GooFit
