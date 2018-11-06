// SourceTimeFunction.cpp
// created by Kuangdai on 7-Apr-2016 
// source time function

#include "SourceTimeFunction.h"

SourceTimeFunction::SourceTimeFunction(const std::vector<Real> &stf, double dt, double shift):
mSTF(stf), mDeltaT(dt), mShift(shift) {
    // nothing
}
