// SourceTimeFunction.h
// created by Kuangdai on 7-Apr-2016 
// source time function


#pragma once

#include <vector>
#include "global.h"

class SourceTimeFunction {
public:
    SourceTimeFunction(const std::vector<Real> &stf, double dt, double shift);
    
    int getSize() const {return mSTF.size();};
    Real getFactor(int tstep) const {return mSTF[tstep];};
    double getDeltaT() const {return mDeltaT;};
    double getShift() const {return mShift;};
    
private:
    std::vector<Real> mSTF;
    double mDeltaT;
    double mShift;
};