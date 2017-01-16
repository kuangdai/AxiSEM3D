// RickerSTF.h
// created by Alex on 8-May-2016
// Ricker stf

#pragma once

#include "STF.h"

class RickerSTF: public STF {
public:
    RickerSTF(double dt, double length, double hdur, double decay);
    std::string verbose() const;

private:
    double mHalfDuration;
    double mDecay;
};
