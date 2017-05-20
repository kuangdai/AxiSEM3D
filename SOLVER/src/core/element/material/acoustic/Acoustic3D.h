// Acoustic3D.h
// created by Kuangdai on 23-Apr-2016 
// 3D acoustic

#pragma once

#include "Acoustic.h"
#include "eigenc.h"

class Acoustic3D: public Acoustic {
public:
    // constructor
    Acoustic3D(const RDMatXN &KFluid);
    
    // STEP 2: strain ==> stress
    void strainToStress(FluidResponse &response) const;
    
    // verbose
    std::string verbose() const {return "Acoustic3D";};
    bool is1D() const {return false;};

    // check compatibility
    void checkCompatibility(int Nr) const;

private:
    RMatXN mKFlat;
};
