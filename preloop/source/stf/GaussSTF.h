// ErfSTF.h
// created by Alex on 9-May-2016
// error function

#pragma once

#include "STF.h"

class GaussSTF: public STF {
public:
    GaussSTF(double dt, double length, double hdur, double decay);
    std::string verbose() const;

private:
    double mHalfDuration;
    double mDecay;
};
