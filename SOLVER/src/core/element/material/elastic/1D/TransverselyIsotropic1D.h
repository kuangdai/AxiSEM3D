// TransverselyIsotropic1D.h
// created by Kuangdai on 22-Apr-2016 
// transversely isotropic 1D material

#pragma once

#include "Elastic1D.h"
#include "eigenc.h"

class TransverselyIsotropic1D: public Elastic1D {
public:
    // constructor
    TransverselyIsotropic1D(const RMatPP &A, const RMatPP &C, const RMatPP &F, 
        const RMatPP &L, const RMatPP &N, Attenuation1D *att):
        Elastic1D(att) , mA(A), mC(C), mF(F), mL(L), mN(N), mN2(two * N) {};
    
    // STEP 2: strain ==>>> stress
    void strainToStress(SolidResponse &response) const;
    
    // verbose
    std::string verbose() const {return "TransverselyIsotropic1D";};
    
    // need TIso
    bool needTIso() const {return true;};
    
private:
    
    // Cijkl scaled by integral factor
    RMatPP mA;
    RMatPP mC;
    RMatPP mF;
    RMatPP mL;
    RMatPP mN;
    RMatPP mN2;
};
