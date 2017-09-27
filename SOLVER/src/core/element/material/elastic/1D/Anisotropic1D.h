// Anisotropic1D.h
// created by Kuangdai on 25-Sep-2017 
// full anisotropic 1D material

#pragma once

#include "Elastic1D.h"
#include "eigenc.h"

class Anisotropic1D: public Elastic1D {
public:
    // constructor
    Anisotropic1D(
        const RMatPP &C11, const RMatPP &C12, const RMatPP &C13, const RMatPP &C14, const RMatPP &C15, const RMatPP &C16,
        const RMatPP &C22, const RMatPP &C23, const RMatPP &C24, const RMatPP &C25, const RMatPP &C26,
        const RMatPP &C33, const RMatPP &C34, const RMatPP &C35, const RMatPP &C36,
        const RMatPP &C44, const RMatPP &C45, const RMatPP &C46,
        const RMatPP &C55, const RMatPP &C56,
        const RMatPP &C66, Attenuation1D *att): Elastic1D(att), 
        mC11(C11), mC12(C12), mC13(C13), mC14(C14), mC15(C15), mC16(C16),
        mC22(C22), mC23(C23), mC24(C24), mC25(C25), mC26(C26), 
        mC33(C33), mC34(C34), mC35(C35), mC36(C36), 
        mC44(C44), mC45(C45), mC46(C46), 
        mC55(C55), mC56(C56),
        mC66(C66) {};
    
    // STEP 2: strain ==>>> stress
    void strainToStress(SolidResponse &response) const;
    
    // verbose
    std::string verbose() const {return "Anisotropic1D";};
    
    // need TIso
    bool needTIso() const {return true;};
    
private:
    
    // Cijkl scaled by integral factor
    RMatPP mC11;
    RMatPP mC12;
    RMatPP mC13;
    RMatPP mC14;
    RMatPP mC15;
    RMatPP mC16;
    
    RMatPP mC22;
    RMatPP mC23;
    RMatPP mC24;
    RMatPP mC25;
    RMatPP mC26;
    
    RMatPP mC33;
    RMatPP mC34;
    RMatPP mC35;
    RMatPP mC36;
    
    RMatPP mC44;
    RMatPP mC45;
    RMatPP mC46;

    RMatPP mC55;
    RMatPP mC56;
    
    RMatPP mC66;
};
