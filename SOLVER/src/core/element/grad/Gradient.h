// Gradient.h
// created by Kuangdai on 19-May-2017 
// elemental gradient

#pragma once

#include "eigenc.h"
#include "eigenp.h"

class FluidElementResponse;
class SolidElementResponse;

class Gradient {
public:
    Gradient(const RDMatPP &dsdxii, const RDMatPP &dsdeta,  
             const RDMatPP &dzdxii, const RDMatPP &dzdeta, 
             const RDMatPP &inv_s, bool axial);
    ~Gradient() {};
    
private:
    
    void computeGrad(FluidElementResponse &response) const;
    void computeQuad(FluidElementResponse &response) const;
    
    void computeGrad9(SolidElementResponse &response) const;
    void computeQuad9(SolidElementResponse &response) const;
    
    void computeGrad6(SolidElementResponse &response) const;
    void computeQuad6(SolidElementResponse &response) const;
    
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