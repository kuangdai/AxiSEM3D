// Attenuation3D_CG4.h
// created by Kuangdai on 29-Apr-2016 
// 3D attenuation on coarse grid 

#pragma once

#include "Attenuation3D.h"
#include "eigen_cg4.h"
#include <vector>

class Attenuation3D_CG4: public Attenuation3D {
public:
    
    Attenuation3D_CG4(int nsls, 
        const RColX &alpha, const RColX &beta, const RColX &gamma, 
        const RMatX4 &dkappa, const RMatX4 &dmu, bool doKappa);
    
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
    RMatX46 mStressR;
    RMatX46 mStrain4;
    std::vector<RMatX46> mMemVar;
    // modules
    RMatX4 mDKappa3; // dkappa * 3
    RMatX4 mDMu;     // dmu
    RMatX4 mDMu2;     // dmu
    bool mDoKappa;
};
