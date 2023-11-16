#pragma once

#include <goofit/PDFs/physics/Amp4BodyGlobals.h>
#include <goofit/PDFs/ParameterContainer.h>

namespace GooFit {

class Cache_SF_LS_TD_EntryFinder final {
  public:
    static const unsigned int _PC_CACHE_IDX   = 5;
    static const unsigned int _PC_TOT_LS_IDX  = 6;
    static const unsigned int _PC_TOT_SF_IDX  = 7;
    static const unsigned int _PC_TOT_AMP_IDX = 8;

    const unsigned int _N_PERM;
    const unsigned int _AMP_IDX;

    Cache_SF_LS_TD_EntryFinder(unsigned int nPerm, unsigned int ampIdx)
        : _N_PERM(nPerm)
        , _AMP_IDX(ampIdx) {}

    void setDalitzId(int idx) { _dalitzFuncId = idx; }

    __device__ inline fpcomplex getLSVal(const ParameterContainer &pc, unsigned int evtNum, unsigned int lsNum) const {
        const unsigned int CACHE_TO_USE = getCacheToUse(pc);
        const unsigned int LS_IND       = getLSInd(pc, evtNum, lsNum);
        return cResSF_TD[CACHE_TO_USE][LS_IND];
    }

    __device__ inline fptype getSFVal(const ParameterContainer &pc, unsigned int evtNum, unsigned int sfNum) const {
        const unsigned int CACHE_TO_USE = getCacheToUse(pc);
        const unsigned int SF_IND       = getSFInd(pc, evtNum, sfNum);
        return (cResSF_TD[CACHE_TO_USE][SF_IND]).real();
    }

  private:
    unsigned int _dalitzFuncId;

    __device__ inline unsigned int
    getSFInd(const ParameterContainer &pc, unsigned int evtNum, unsigned int sfNum) const {
        const unsigned int NUM_LS    = getNumLS(pc);
        const unsigned int OFFSET    = getOffset(pc);
        const unsigned int TOTAL_LS  = getTotLS(pc);
        const unsigned int TOTAL_AMP = getTotAMP(pc);
        unsigned int idx = AmpIndices[TOTAL_AMP + _AMP_IDX + 4 + NUM_LS + sfNum];
        return evtNum * OFFSET + TOTAL_LS + idx;
    }

    __device__ inline unsigned int
    getLSInd(const ParameterContainer &pc, unsigned int evtNum, unsigned int lsNum) const {
        const unsigned int TOTAL_AMP = getTotAMP(pc);
        const unsigned int OFFSET    = getOffset(pc);
        unsigned int idx = AmpIndices[TOTAL_AMP + _AMP_IDX + 4 + lsNum];
        return evtNum * OFFSET + idx;
    }

    __device__ static inline unsigned int getCacheToUse(const ParameterContainer &pc) { return pc.getConstant(_PC_CACHE_IDX); }

    __device__ inline unsigned int getSFStep(const ParameterContainer &pc) const { return getNumSF(pc) / _N_PERM; }

    __device__ inline unsigned int getLSStep(const ParameterContainer &pc) const { return getNumLS(pc) / _N_PERM; }

    __device__ static inline unsigned int getTotLS(const ParameterContainer &pc) { return pc.getConstant(_PC_TOT_LS_IDX); }

    __device__ static inline unsigned int getTotSF(const ParameterContainer &pc) { return pc.getConstant(_PC_TOT_SF_IDX); }

    __device__ static inline unsigned int getTotAMP(const ParameterContainer &pc) { return pc.getConstant(_PC_TOT_AMP_IDX); }

    __device__ inline unsigned int getNumSF(const ParameterContainer &pc) const {
        const unsigned int totalAMP = getTotAMP(pc);
        return AmpIndices[totalAMP + _AMP_IDX + 1];
    }

    __device__ inline unsigned int getNumLS(const ParameterContainer &pc) const {
        const unsigned int totalAMP = getTotAMP(pc);
        return AmpIndices[totalAMP + _AMP_IDX];
    }

    __device__ static inline unsigned int getOffset(const ParameterContainer &pc) { return getTotLS(pc) + getTotSF(pc); }

    friend class AmpCalc_TD;
    friend struct SF_EntryDebugger_TD;
    friend struct LS_EntryDebugger_TD;
};

} // namespace GooFit