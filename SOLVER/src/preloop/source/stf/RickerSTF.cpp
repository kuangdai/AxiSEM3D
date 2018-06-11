// RickerSTF.cpp
// created by Alex on 8-May-2016
// Ricker stf

#include "RickerSTF.h"
#include <cmath>
#include <string>
#include <sstream>

RickerSTF::RickerSTF(double dt, double duration, double hdur, double decay):
mHalfDuration(hdur), mDecay(decay) {
    mDeltaT = dt;
    int nStepBeforeZero = ceil(2.5 * mHalfDuration / mDeltaT);
    int nStepAfterZero = ceil(duration / mDeltaT);
    mShift = nStepBeforeZero * mDeltaT;
    int nStep = nStepBeforeZero + nStepAfterZero;
    for (int i = 0; i <= nStep; i++) {
        double t = -mShift + i * mDeltaT;
        mSTF.push_back(2. * pow(mDecay / mHalfDuration, 2.) * exp(-pow((mDecay / mHalfDuration * t), 2.)) 
            * (2. * pow(t, 2.) * pow(mDecay / mHalfDuration, 2.) - 1.));
    }
}

std::string RickerSTF::verbose() const {
    std::stringstream ss;
    ss << "\n=================== Source Time Function ===================" << std::endl;
    ss << "  Time Step               =   " << mDeltaT << std::endl;
    ss << "  Number of Steps         =   " << mSTF.size() << std::endl;
    ss << "  Total Duration          =   " << mDeltaT * mSTF.size() << std::endl;
    ss << "  Duration after Origin   =   " << mDeltaT * mSTF.size() - mShift << std::endl;
    ss << "  Shift before Origin     =   " << mShift << std::endl;
    ss << "  Time Series Type        =   Ricker" << std::endl;
    ss << "  Half Duration           =   " << mHalfDuration << std::endl;
    ss << "  Decay Factor            =   " << mDecay << std::endl;
    ss << "=================== Source Time Function ===================\n" << std::endl;
    return ss.str();
}
