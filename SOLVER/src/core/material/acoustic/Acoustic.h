// Acoustic.h
// created by Kuangdai on 23-Apr-2016 
// base class of acoustic constitutive relation

#pragma once

#include "eigenc.h"

class Acoustic {
public:
    virtual ~Acoustic() {};
    
    // STEP 2: strain ==> stress
    virtual void strainToStress(const vec_ar3_CMatPP &strain, vec_ar3_CMatPP &stress, int Nu) const = 0;
        
    // check compatibility
    virtual void checkCompatibility(int Nr) const {}; 
    
    // verbose
    virtual std::string verbose() const = 0;
    
};

