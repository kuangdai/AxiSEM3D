// Elastic3D.h
// created by Kuangdai on 29-Apr-2016 
// base class of 3D elasticity

#pragma once

#include "Elastic.h"

class Attenuation3D;

class Elastic3D: public Elastic {
public:
    Elastic3D(Attenuation3D *att);
    virtual ~Elastic3D();
    
    // check compatibility
    virtual void checkCompatibility(int Nr) const; 
    
    // reset to zero 
    void resetZero(); 
    
    // 1D or Fourier space
    bool is1D() const {return false;};
    
protected:
    Attenuation3D *mAttenuation;
};
