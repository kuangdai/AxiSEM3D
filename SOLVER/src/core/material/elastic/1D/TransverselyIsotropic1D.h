// TransverselyIsotropic1D.h
// created by Kuangdai on 22-Apr-2016 
// transversely isotropic 1D material

#pragma once

#include "Elastic1D.h"
#include "eigenp.h"

class TransverselyIsotropic1D: public Elastic1D {
public:
    // constructor
    TransverselyIsotropic1D(const RDMatPP &theta, const RMatPP &A, const RMatPP &C, 
        const RMatPP &F, const RMatPP &L, const RMatPP &N, Attenuation1D *att);
    
    // STEP 2: strain ==>>> stress
    void strainToStress(const vec_ar9_CMatPP &strain, vec_ar9_CMatPP &stress, int Nu) const;
    
    // verbose
    std::string verbose() const {return "TransverselyIsotropic1D";};
    
private:    
    void rotateStrainToTIso(const ar9_CMatPP &strainCyln, ar6_CMatPP &strainTIso) const;
    void rotateStressToCyln(const ar6_CMatPP &stressTIso, ar9_CMatPP &stressCyln) const;

private:
    
    // orientation of axis of symmetry of transverse isotropy
    RMatPP mSin1t = RMatPP::Zero();
    RMatPP mCos1t = RMatPP::Zero();
    RMatPP mSin2t = RMatPP::Zero();
    RMatPP mCos2t = RMatPP::Zero();
    
    // Cijkl scaled by integral factor
    RMatPP mA;
    RMatPP mC;
    RMatPP mF;
    RMatPP mL;
    RMatPP mN;
    RMatPP mN2;
};
