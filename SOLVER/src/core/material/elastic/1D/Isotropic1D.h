// Isotropic1D.h
// created by Kuangdai on 3-Apr-2016 
// isotropic 1D material

#pragma once

#include "Elastic1D.h"

class Isotropic1D: public Elastic1D {
public:
    // constructor
    Isotropic1D(const RMatPP &lambda, const RMatPP &mu, Attenuation1D *att);
    
    // STEP 2: strain ==>>> stress
    void strainToStress(const vec_ar9_CMatPP &strain, vec_ar9_CMatPP &stress, int Nu) const;
    
    // verbose
    std::string verbose() const {return "Isotropic1D";};
                    
private:
    // Cijkl scaled by integral factor
    RMatPP mLambda; 
    RMatPP mMu; 
    RMatPP mMu2;
    
};
