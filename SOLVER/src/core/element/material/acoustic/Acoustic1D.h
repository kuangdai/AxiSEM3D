// Acoustic1D.h
// created by Kuangdai on 23-Apr-2016 
// 1D acoustic

#pragma once

#include "Acoustic.h"
#include "eigenc.h"

class Acoustic1D: public Acoustic {
public:
    // constructor
    Acoustic1D(const RDMatXN &KFluid);
    
    // STEP 2: strain ==> stress
    void strainToStress(FluidElementResponse &response) const;
    
    // verbose
    std::string verbose() const {return "Acoustic1D";};
    bool is1D() const {return true;};
    
private:
    RMatPP mKStruct; 
};
