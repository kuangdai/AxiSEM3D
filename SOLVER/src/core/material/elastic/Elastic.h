// Elastic.h
// created by Kuangdai on 2-Apr-2016 
// base class of elasticity

#pragma once

#include "eigenc.h"

class Elastic {
public:
    
    virtual ~Elastic() {};
    
    // STEP 2: strain ==> stress
    virtual void strainToStress(const vec_ar9_CMatPP &strain, vec_ar9_CMatPP &stress, int Nu) const = 0;
        
    // check compatibility
    virtual void checkCompatibility(int Nr, bool isVoigt) const = 0; 
    
    // verbose
    virtual std::string verbose() const = 0;
    
    // reset to zero 
    virtual void resetZero() = 0; 
    
};
