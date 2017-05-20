// Gradient.h
// created by Kuangdai on 19-May-2017 
// elemental gradient

#pragma once

#include "eigenc.h"
#include "eigenp.h"

class Gradient {
public:
    Gradient(const RDMatPP &dsdxii, const RDMatPP &dsdeta,  
             const RDMatPP &dzdxii, const RDMatPP &dzdeta, 
             const RDMatPP &inv_s, bool axial);
    ~Gradient() {};
    
    void computeGrad(const vec_CMatPP &u, vec_ar3_CMatPP &u_i, int Nu, int nyquist) const;
    void computeQuad(vec_CMatPP &f, const vec_ar3_CMatPP &f_i, int Nu, int nyquist) const;
    
    void computeGrad9(const vec_ar3_CMatPP &ui, vec_ar9_CMatPP &ui_j, int Nu, int nyquist) const;
    void computeQuad9(vec_ar3_CMatPP &fi, const vec_ar9_CMatPP &fi_j, int Nu, int nyquist) const;
    
    void computeGrad6(const vec_ar3_CMatPP &ui, vec_ar6_CMatPP &eij, int Nu, int nyquist) const;
    void computeQuad6(vec_ar3_CMatPP &fi, const vec_ar6_CMatPP &sij, int Nu, int nyquist) const;

private:    
    // operators
    RMatPP mDsDxii;
    RMatPP mDsDeta;
    RMatPP mDzDxii;
    RMatPP mDzDeta;
    RMatPP mInv_s;
    
    // axis
    bool mAxial;
    RMatPP *sG_xii;
    RMatPP *sGT_xii;
    RMatPP *sG_eta;
    RMatPP *sGT_eta;
    
//-------------------------- static --------------------------//
public: 
    // set G Mat, shared by all elements
    static void setGMat(const RDMatPP &G_GLL, const RDMatPP &G_GLJ);
    
private: 
    // G Mat storage
    static RMatPP sG_GLL;
    static RMatPP sG_GLJ;
    static RMatPP sGT_GLL;
    static RMatPP sGT_GLJ;
};