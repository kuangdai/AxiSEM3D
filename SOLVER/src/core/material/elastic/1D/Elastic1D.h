// Elastic1D.h
// created by Kuangdai on 29-Apr-2016 
// base class of 1D elasticity

#pragma once

#include "eigenc.h"
#include "Elastic.h"

class Attenuation1D;

class Elastic1D: public Elastic {
public:
    Elastic1D(Attenuation1D *att);
    virtual ~Elastic1D();
    
    // check compatibility
    virtual void checkCompatibility(int Nr, bool isVoigt) const; 
    
    // reset to zero 
    void resetZero(); 
    
protected:
    Attenuation1D *mAttenuation;
};
