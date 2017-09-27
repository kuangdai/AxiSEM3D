// Anisotropic3D.h
// created by Kuangdai on 25-Sep-2017 
// full anisotropic 3D material 

#pragma once

#include "Elastic3D.h"
#include "eigenc.h"

class Anisotropic3D: public Elastic3D {
public:
    // constructor
    Anisotropic3D(
        const RMatXN &C11, const RMatXN &C12, const RMatXN &C13, const RMatXN &C14, const RMatXN &C15, const RMatXN &C16,
        const RMatXN &C22, const RMatXN &C23, const RMatXN &C24, const RMatXN &C25, const RMatXN &C26,
        const RMatXN &C33, const RMatXN &C34, const RMatXN &C35, const RMatXN &C36,
        const RMatXN &C44, const RMatXN &C45, const RMatXN &C46,
        const RMatXN &C55, const RMatXN &C56,
        const RMatXN &C66, Attenuation3D *att): Elastic3D(att), 
        mC11(C11), mC12(C12), mC13(C13), mC14(C14), mC15(C15), mC16(C16),
        mC22(C22), mC23(C23), mC24(C24), mC25(C25), mC26(C26), 
        mC33(C33), mC34(C34), mC35(C35), mC36(C36), 
        mC44(C44), mC45(C45), mC46(C46), 
        mC55(C55), mC56(C56),
        mC66(C66) {};
        
    // STEP 2: strain ==>>> stress
    void strainToStress(SolidResponse &response) const;
    
    // check compatibility
    void checkCompatibility(int Nr) const; 
    
    // verbose
    std::string verbose() const {return "Anisotropic3D";};
    
    // need TIso
    bool needTIso() const {return true;};
    
private:
    // Cijkl scaled by integral factor
    RMatXN mC11;
    RMatXN mC12;
    RMatXN mC13;
    RMatXN mC14;
    RMatXN mC15;
    RMatXN mC16;
    
    RMatXN mC22;
    RMatXN mC23;
    RMatXN mC24;
    RMatXN mC25;
    RMatXN mC26;
    
    RMatXN mC33;
    RMatXN mC34;
    RMatXN mC35;
    RMatXN mC36;
    
    RMatXN mC44;
    RMatXN mC45;
    RMatXN mC46;

    RMatXN mC55;
    RMatXN mC56;
    
    RMatXN mC66;
};
