// TransverselyIsotropic3D.h
// created by Kuangdai on 22-Apr-2016 
// transversely isotropic 3D material 

#pragma once

#include "Elastic3D.h"
#include "eigenc.h"

class TransverselyIsotropic3D: public Elastic3D {
public:
    // constructor
    TransverselyIsotropic3D(const RMatXN &A, const RMatXN &C, const RMatXN &F, 
        const RMatXN &L, const RMatXN &N, Attenuation3D *att):
        Elastic3D(att) , mA(A), mC(C), mF(F), mL(L), mN(N), mN2(two * N) {};
        
    // STEP 2: strain ==>>> stress
    void strainToStress(SolidResponse &response) const;
    
    // check compatibility
    void checkCompatibility(int Nr) const; 
    
    // verbose
    std::string verbose() const {return "TransverselyIsotropic3D";};
    
    // need TIso
    bool needTIso() const {return true;};
    
private:
    // Cijkl scaled by integral factor
    RMatXN mA;
    RMatXN mC;
    RMatXN mF;
    RMatXN mL;
    RMatXN mN;
    RMatXN mN2;
};
