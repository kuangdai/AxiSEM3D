// ErfSTF.h
// created by Kuangdai on 9-May-2016 
// error function

#pragma once

#include "STF.h"

class ErfSTF: public STF {
public: 
    ErfSTF(double dt, double length, double hdur, double decay);
    std::string verbose() const;

private:    
    double mHalfDuration;
    double mDecay;
};