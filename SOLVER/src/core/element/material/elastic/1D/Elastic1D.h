// Elastic1D.h
// created by Kuangdai on 29-Apr-2016 
// base class of 1D elasticity

#pragma once

#include "Elastic.h"

class Attenuation1D;

class Elastic1D: public Elastic {
public:
    Elastic1D(Attenuation1D *att);
    virtual ~Elastic1D();
    
    // check compatibility
    virtual void checkCompatibility(int Nr) const; 
    
    // reset to zero 
    void resetZero(); 
    
    // 1D or Fourier space
    bool is1D() const {return true;};
    
protected:
    Attenuation1D *mAttenuation;
};
