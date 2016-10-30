// MassOcean3D.h
// created by Kuangdai on 2-Jun-2016 
// 3D mass with ocean

#pragma once

#include "Mass.h"
#include "eigenp.h"

class MassOcean3D : public Mass {
public:
    MassOcean3D(const RDColX &mass, const RDColX &massOcean, const RDMatX3 &normal);
    
    // compute accel in-place
    void computeAccel(CMatX3 &stiff) const;
    void computeAccel(CColX &stiff) const;
    
    void checkCompatibility(int nr) const;
    
    // verbose
    std::string verbose() const {return "MassOcean3D";};
    
private:
    RColX mInvMass;
    RMatX3 mNormal_scal; 
};
