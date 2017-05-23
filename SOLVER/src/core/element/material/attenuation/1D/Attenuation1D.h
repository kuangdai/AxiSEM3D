// Attenuation1D.h
// created by Kuangdai on 29-Apr-2016 
// base class of 1D attenuation 

#pragma once

#include "Attenuation.h"
#include "eigenp.h"

class Attenuation1D: public Attenuation {
public:
    
    Attenuation1D(int nsls, const RColX &alpha, const RColX &beta, const RColX &gamma):
        Attenuation(nsls, alpha, beta, gamma) {};
    virtual ~Attenuation1D() {};
        
    // STEP 2.1: R ==> stress
    virtual void applyToStress(vec_ar6_CMatPP &stress) const = 0;

    // STEP 2.3: strain ==> R
    virtual void updateMemoryVariables(const vec_ar6_CMatPP &strain) = 0;
};
