// PreloopGradient.h
// created by Kuangdai on 20-Sep-2016 
// elemental gradient

#pragma once

#include "eigenp.h"

class PreloopGradient {
public:    
    
    PreloopGradient(const RDMatPP &dsdxii, const RDMatPP &dsdeta,  
                    const RDMatPP &dzdxii, const RDMatPP &dzdeta, 
                    const RDMatPP &inv_s, bool axial);
    void gradScalar(const vec_CDMatPP &u, vec_ar3_CDMatPP &u_i, int Nu, int nyquist) const;
    
private:
    // axial
    bool mAxial;
    
    // operators
    RDMatPP mDsDxii;
    RDMatPP mDsDeta;
    RDMatPP mDzDxii;
    RDMatPP mDzDeta;
    RDMatPP mInv_s;
    
//-------------------------- static --------------------------//
public: 
    // set G Mat, shared by all elements
    static void setGMat(const RDMatPP &G_GLL, const RDMatPP &G_GLJ);
    
protected: 
    // G Mat storage
    static RDMatPP sG_GLL;
    static RDMatPP sG_GLJ;
    static RDMatPP sGT_GLL;
    static RDMatPP sGT_GLJ;
};