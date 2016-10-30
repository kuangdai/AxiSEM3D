// Mass1D.h
// created by Kuangdai on 3-Apr-2016 
// 1D mass

#pragma once

#include "Mass.h"

class Mass1D : public Mass {
public:
    Mass1D(Real invMass);

    // compute accel in-place
    void computeAccel(CMatX3 &stiff) const;
    void computeAccel(CColX &stiff) const;
    
    // verbose
    std::string verbose() const {return "Mass1D";};
    
private:
    // the scalar mass
    Real mInvMass;
};
