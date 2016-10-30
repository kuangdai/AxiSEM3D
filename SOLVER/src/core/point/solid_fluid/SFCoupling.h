// SFCoupling.h
// created by Kuangdai on 5-Apr-2016 
// solid-fluid boundary condition

#pragma once

#include "eigenc.h"

class SFCoupling {
public:   
    virtual ~SFCoupling() {};
    
    // solid-fluid coupling
    virtual void coupleFluidToSolid(const CColX &fluidStiff, CMatX3 &solidStiff) const = 0; 
    virtual void coupleSolidToFluid(const CMatX3 &solidDispl, CColX &fluidStiff) const = 0;
    
    // verbose
    virtual std::string verbose() const = 0;    
    
    // check compatibility
    virtual void checkCompatibility(int nr) const {};
};
