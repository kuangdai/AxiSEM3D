// ErfSTF.cpp
// created by Kuangdai on 9-May-2016 
// error function

#include "ErfSTF.h"
#include <cmath>
#include <string>
#include <sstream>

ErfSTF::ErfSTF(double dt, double duration, double hdur, double decay): 
mHalfDuration(hdur), mDecay(decay) {
    mDeltaT = dt;
    int nStepBeforeZero = ceil(2.5 * mHalfDuration / mDeltaT);
    int nStepAfterZero = ceil(duration / mDeltaT);
    mShift = nStepBeforeZero * mDeltaT;
    int nStep = nStepBeforeZero + nStepAfterZero;
    for (int i = 0; i <= nStep; i++) {
        double t = -mShift + i * mDeltaT;
        mSTF.push_back(erf(mDecay / mHalfDuration * t) * 0.5 + 0.5);
    }
}

std::string ErfSTF::verbose() const {
    std::stringstream ss;
    ss << "\n=================== Source Time Function ===================" << std::endl;
    ss << "  Time Step               =   " << mDeltaT << std::endl;
    ss << "  Number of Steps         =   " << mSTF.size() << std::endl;
    ss << "  Total Duration          =   " << mDeltaT * mSTF.size() << std::endl;
    ss << "  Duration after Origin   =   " << mDeltaT * mSTF.size() - mShift << std::endl;
    ss << "  Shift before Origin     =   " << mShift << std::endl;
    ss << "  Time Series Type        =   Erf" << std::endl;
    ss << "  Half Duration           =   " << mHalfDuration << std::endl;
    ss << "  Decay Factor            =   " << mDecay << std::endl;
    ss << "=================== Source Time Function ===================\n" << std::endl;
    return ss.str();
}
