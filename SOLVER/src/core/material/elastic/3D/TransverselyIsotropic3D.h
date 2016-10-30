// TransverselyIsotropic3D.h
// created by Kuangdai on 22-Apr-2016 
// transversely isotropic 3D material 

#pragma once

#include "Elastic3D.h"
#include "eigenp.h"

class TransverselyIsotropic3D: public Elastic3D {
public:

    TransverselyIsotropic3D(const RDRowN &theta, 
        const RMatXN &A, const RMatXN &C, const RMatXN &F, 
        const RMatXN &L, const RMatXN &N, Attenuation3D *att);
        
    // STEP 2: strain ==>>> stress
    void strainToStress(const vec_ar9_CMatPP &strain, vec_ar9_CMatPP &stress, int Nu) const;
    
    // check compatibility
    void checkCompatibility(int Nr, bool isVoigt) const; 
    
    // verbose
    std::string verbose() const {return "TransverselyIsotropic3D";};
    
private:    
    
    void rotateStrainToTIso(const CMatXN6 &strainCyln, CMatXN6 &strainTIso, int Nu) const;
    void rotateStressToCyln(const CMatXN6 &stressTIso, CMatXN6 &stressCyln, int Nu) const;
                
private:
    
    // orientation of axis of symmetry of transverse isotropy
    RMatXN mSin1t;
    RMatXN mCos1t;
    RMatXN mSin2t;
    RMatXN mCos2t;
    
    // Cijkl scaled by integral factor
    RMatXN mA;
    RMatXN mC;
    RMatXN mF;
    RMatXN mL;
    RMatXN mN;
    RMatXN mN2;
};
