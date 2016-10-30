// Isotropic3D.h
// created by Kuangdai on 22-Apr-2016 
// isotropic 3D material

#pragma once

#include "Elastic3D.h"

class Isotropic3D: public Elastic3D {
public:
    
    Isotropic3D(const RMatXN &lambda, const RMatXN &mu, Attenuation3D *att);
    
    // STEP 2: strain ==>>> stress
    void strainToStress(const vec_ar9_CMatPP &strain, vec_ar9_CMatPP &stress, int Nu) const;
    
    // check compatibility
    void checkCompatibility(int Nr, bool isVoigt) const; 
    
    // verbose
    std::string verbose() const {return "Isotropic3D";};
                
private:    
    // Cijkl scaled by integral factor
    RMatXN mLambda;
    RMatXN mMu;
    RMatXN mMu2;
};
