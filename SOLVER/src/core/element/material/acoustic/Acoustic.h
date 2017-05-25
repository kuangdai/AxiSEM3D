// Acoustic.h
// created by Kuangdai on 23-Apr-2016 
// base class of acoustic constitutive relation

#pragma once

class FluidResponse;
#include <string>
#include "eigenp.h"

class Acoustic {
public:
    virtual ~Acoustic() {};
    
    // STEP 2: strain ==> stress
    virtual void strainToStress(FluidResponse &response) const = 0;
    
    // verbose
    virtual std::string verbose() const = 0;
    
    // 1D or Fourier space
    virtual bool is1D() const = 0;

    // check compatibility
    virtual void checkCompatibility(int Nr) const {}; 
};

