// Attenuation1D_Full.h
// created by Kuangdai on 29-Apr-2016 
// 1D attenuation on full grid

#pragma once

#include "Attenuation1D.h"
#include <vector>

class Attenuation1D_Full: public Attenuation1D {
public:
    
    Attenuation1D_Full(int nsls, const RColX &alpha, 
        const RColX &beta, const RColX &gamma, int Nu, 
        const RMatPP &dkappa, const RMatPP &dmu, bool doKappa);
        
    // STEP 2.1: R ==> stress
    void applyToStress(vec_ar6_CMatPP &stress) const;

    // STEP 2.3: strain ==> R
    void updateMemoryVariables(const vec_ar6_CMatPP &strain);
    
    // check memory variable size
    void checkCompatibility(int Nr) const;
    
    // reset to zero 
    void resetZero(); 
    
private:
    // memory variables
    vec_ar6_CMatPP mStressR;
    std::vector<vec_ar6_CMatPP> mMemVar;
    // modules
    RMatPP mDKappa3; // dkappa * 3
    RMatPP mDMu;     // dmu
    RMatPP mDMu2;    // dmu * 2
    bool mDoKappa;
};

