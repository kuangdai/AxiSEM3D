// GradientVoigt.h
// created by Kuangdai on 8-May-2016 
// elemental gradient using Voigt notation

#pragma once

#include "Gradient.h"

class GradientVoigt: public Gradient {
public:    
    
    GradientVoigt(const RMatPP &dsdxii, const RMatPP &dsdeta,  
        const RMatPP &dzdxii, const RMatPP &dzdeta, const RMatPP &inv_s);
        
    // Voigt notation, strain 6
    void gradVector(const vec_ar3_CMatPP &ui, vec_ar9_CMatPP &eij, int Nu, int nyquist) const;
    void quadVector(const vec_ar9_CMatPP &sij, vec_ar3_CMatPP &fi, int Nu, int nyquist) const;
    
    // if Voigt notation is used
    bool isVoigt() const {return true;};
};

