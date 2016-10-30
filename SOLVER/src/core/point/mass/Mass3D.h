// Mass3D.h
// created by Kuangdai on 2-Jun-2016 
// 3D mass

#pragma once

#include "Mass.h"

class Mass3D : public Mass {
public:
    Mass3D(const RColX &invMass);
    
    // compute accel in-place
    void computeAccel(CMatX3 &stiff) const;
    void computeAccel(CColX &stiff) const;
    
    void checkCompatibility(int nr) const;
    
    // verbose
    std::string verbose() const {return "Mass3D";};
    
private:
    RColX mInvMass;
};
