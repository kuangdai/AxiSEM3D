// Isotropic3D.h
// created by Kuangdai on 22-Apr-2016 
// isotropic 3D material

#pragma once

#include "Elastic3D.h"
#include "eigenc.h"

class Isotropic3D: public Elastic3D {
public:
    // constructor
    Isotropic3D(const RMatXN &lambda, const RMatXN &mu, Attenuation3D *att):
        Elastic3D(att), mLambda(lambda), mMu(mu), mMu2(two * mu) {};
    
    // STEP 2: strain ==>>> stress
    void strainToStress(SolidResponse &response) const;
    
    // check compatibility
    void checkCompatibility(int Nr) const; 
    
    // verbose
    std::string verbose() const {return "Isotropic3D";};
    
    // need TIso
    bool needTIso() const {return false;};
                
private:    
    // Cijkl scaled by integral factor
    RMatXN mLambda;
    RMatXN mMu;
    RMatXN mMu2;
};
