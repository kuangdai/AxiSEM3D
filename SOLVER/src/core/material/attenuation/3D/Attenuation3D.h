// Attenuation3D.h
// created by Kuangdai on 29-Apr-2016 
// base class of 3D attenuation 

#pragma once

#include "Attenuation.h"

class Attenuation3D: public Attenuation {
public:
    
    Attenuation3D(int nsls, const RColX &alpha, 
        const RColX &beta, const RColX &gamma);
    virtual ~Attenuation3D() {};
        
    // STEP 2.1: R ==> stress
    virtual void applyToStress(RMatXN6 &stress) const = 0;

    // STEP 2.3: strain ==> R
    virtual void updateMemoryVariables(const RMatXN6 &strain) = 0;
    
    // reset to zero 
    virtual void resetZero() = 0; 
};
