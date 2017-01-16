// GaussSTF.h
// created by Alex on 8-May-2016
// gaussian stf

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
