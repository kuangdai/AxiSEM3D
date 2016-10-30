// Gradient.h
// created by Kuangdai on 23-Apr-2016 
// elemental gradient

#pragma once

#include "eigenc.h"

class Gradient {
public:    
    
    Gradient(const RMatPP &dsdxii, const RMatPP &dsdeta,  
             const RMatPP &dzdxii, const RMatPP &dzdeta, const RMatPP &inv_s);
             
    virtual ~Gradient() {};         
    
    // scalar
    virtual void gradScalar(const vec_CMatPP &u, vec_ar3_CMatPP &u_i, int Nu, int nyquist) const;
    virtual void quadScalar(const vec_ar3_CMatPP &f_i, vec_CMatPP &f, int Nu, int nyquist) const;
    
    // deformation gradient 3x3
    virtual void gradVector(const vec_ar3_CMatPP &ui, vec_ar9_CMatPP &ui_j, int Nu, int nyquist) const;
    virtual void quadVector(const vec_ar9_CMatPP &fi_j, vec_ar3_CMatPP &fi, int Nu, int nyquist) const;
    
    // if Voigt notation is used
    virtual bool isVoigt() const {return false;};
    
protected:
    
    // operators
    RMatPP mDsDxii;
    RMatPP mDsDeta;
    RMatPP mDzDxii;
    RMatPP mDzDeta;
    RMatPP mInv_s;
    
//-------------------------- static --------------------------//
public: 
    // set G Mat, shared by all elements
    static void setGMat(const RMatPP &G_GLL, const RMatPP &G_GLJ);
    
protected: 
    // G Mat storage
    static RMatPP sG_GLL;
    static RMatPP sG_GLJ;
    static RMatPP sGT_GLL;
    static RMatPP sGT_GLJ;
};