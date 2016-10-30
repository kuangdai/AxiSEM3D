// SourceTimeFunction.h
// created by Kuangdai on 7-Apr-2016 
// source time function


#pragma once

#include <vector>
#include "global.h"

class SourceTimeFunction {
public:
    SourceTimeFunction(const std::vector<Real> &stf, Real dt, Real shift);
    
    int getSize() const {return mSTF.size();};
    Real getFactor(int tstep) const {return mSTF[tstep];};
    Real getDeltaT() const {return mDeltaT;};
    Real getShift() const {return mShift;};
    
private:
    std::vector<Real> mSTF;
    Real mDeltaT;
    Real mShift;
};