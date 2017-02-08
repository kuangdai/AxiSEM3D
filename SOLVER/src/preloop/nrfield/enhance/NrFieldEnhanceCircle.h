// NrFieldEnhanceCircle.h
// created by Kuangdai on 2-Feb-2017
// circle-shaped enhanced nr integer field

#pragma once
#include "eigenp.h"
#include "NrFieldEnhance.h"

class Parameters;

class NrFieldEnhanceCircle: public NrFieldEnhance {
public:
    NrFieldEnhanceCircle(int ref, bool decrease, double r, double theta,
        double diameter, double hwhm, double value);
    std::string verbose() const;
    double getValue(const RDCol2 &sz_target) const;
    
protected:
    double mR;
    double mTheta;
    double mDiameter;
    double mHWHM;
    double mValue;
};

