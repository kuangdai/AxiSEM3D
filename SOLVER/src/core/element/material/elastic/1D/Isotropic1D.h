// Isotropic1D.h
// created by Kuangdai on 3-Apr-2016 
// isotropic 1D material

#pragma once

#include "Elastic1D.h"
#include "eigenc.h"

class Isotropic1D: public Elastic1D {
public:
    // constructor
    Isotropic1D(const RMatPP &lambda, const RMatPP &mu, Attenuation1D *att):
        Elastic1D(att), mLambda(lambda), mMu(mu), mMu2(two * mu) {};
    
    // STEP 2: strain ==>>> stress
    void strainToStress(SolidResponse &response) const;
    
    // verbose
    std::string verbose() const {return "Isotropic1D";};
    
    // need TIso
    bool needTIso() const {return false;};
                    
private:
    // Cijkl scaled by integral factor
    RMatPP mLambda; 
    RMatPP mMu; 
    RMatPP mMu2;
};
