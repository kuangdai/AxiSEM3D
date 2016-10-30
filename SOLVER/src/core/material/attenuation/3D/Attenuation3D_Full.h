// Attenuation3D_Full.h
// created by Kuangdai on 29-Apr-2016 
// 3D attenuation on full grid 

#pragma once

#include "Attenuation3D.h"
#include <vector>

class Attenuation3D_Full: public Attenuation3D {
    
public:
    Attenuation3D_Full(int nsls, 
        const RColX &alpha, const RColX &beta, const RColX &gamma, 
        const RMatXN &dkappa, const RMatXN &dmu, bool doKappa);
        
    // STEP 2.1: R ==> stress
    void applyToStress(RMatXN6 &stress) const;

    // STEP 2.3: strain ==> R
    void updateMemoryVariables(const RMatXN6 &strain);
    
    // check memory variable size
    void checkCompatibility(int Nr) const;
    
    // reset to zero 
    void resetZero(); 
    
private:
    // memory variables
    RMatXN6 mStressR;
    std::vector<RMatXN6> mMemVar;
    // modules
    RMatXN mDKappa3; // dkappa * 3
    RMatXN mDMu;     // dmu
    RMatXN mDMu2;    // 2 dmu
    bool mDoKappa;
};
