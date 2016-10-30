// GradientAxial.h
// created by Kuangdai on 23-Apr-2016 
// elemental gradient

#pragma once

#include "Gradient.h"

class GradientAxial: public Gradient {
public:    
    
    GradientAxial(const RMatPP &dsdxii, const RMatPP &dsdeta, 
                  const RMatPP &dzdxii, const RMatPP &dzdeta, const RMatPP &inv_s);
    
    virtual ~GradientAxial() {};
    
    // scalar
    void gradScalar(const vec_CMatPP &u, vec_ar3_CMatPP &u_i, int Nu, int nyquist) const;
    void quadScalar(const vec_ar3_CMatPP &f_i, vec_CMatPP &f, int Nu, int nyquist) const;
    
    // deformation gradient 3x3
    virtual void gradVector(const vec_ar3_CMatPP &ui, vec_ar9_CMatPP &ui_j, int Nu, int nyquist) const;
    virtual void quadVector(const vec_ar9_CMatPP &fi_j, vec_ar3_CMatPP &fi, int Nu, int nyquist) const;
    
    // if Voigt notation is used
    virtual bool isVoigt() const {return false;};
};

