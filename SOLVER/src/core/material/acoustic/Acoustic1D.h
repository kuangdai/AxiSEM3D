// Acoustic1D.h
// created by Kuangdai on 23-Apr-2016 
// 1D acoustic

#pragma once

#include "Acoustic.h"

class Acoustic1D: public Acoustic {
public:
    // constructor
    Acoustic1D(const RMatPP &KFluid);
    
    // STEP 2: strain ==> stress
    void strainToStress(const vec_ar3_CMatPP &strain, vec_ar3_CMatPP &stress, int Nu) const;
    
    // verbose
    std::string verbose() const {return "Acoustic1D";};

private:
    
    RMatPP mK; 

};