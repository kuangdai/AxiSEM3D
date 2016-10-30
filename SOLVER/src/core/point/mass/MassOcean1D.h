// MassOcean1D.h
// created by Kuangdai on 3-Apr-2016 
// 1D mass with ocean

#pragma once

#include "Mass.h"

class MassOcean1D : public Mass {
public:
    MassOcean1D(double mass, double massOcean, double theta);

    // compute accel in-place
    void computeAccel(CMatX3 &stiff) const;
    void computeAccel(CColX &stiff) const;
    
    // verbose
    std::string verbose() const {return "MassOcean1D";};
    
private:
    // the scalar mass
    Real mInvMassZ;
    Real mInvMassR;
    Real mSint;
    Real mCost;
};
