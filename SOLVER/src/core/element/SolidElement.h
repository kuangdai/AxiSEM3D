// SolidElement.h
// created by Kuangdai on 29-Mar-2016 
// base class of solid element

#pragma once

#include "Element.h"

class Elastic;
class CrdTransTIsoSolid;

// static workspaces
struct SolidResponse {
    // disp
    vec_ar3_CMatPP mDispl;
    // strain
    vec_ar6_CMatPP mStrain6;
    vec_ar9_CMatPP mStrain9;
    // stress
    vec_ar6_CMatPP mStress6;
    vec_ar9_CMatPP mStress9;
    // stiff
    vec_ar3_CMatPP mStiff;
    // size
    int mNu = 0;
    int mNr = 1;
    int mNyquist = 0;
    void setNr(int nr) {
        mNr = nr;
        mNu = nr / 2;
        mNyquist = (int)(mNr % 2 == 0);
    };
};

class SolidElement : public Element {
    friend class FluidElement;
    
public:
    
    SolidElement(Gradient *grad, PRT *prt, const std::array<Point *, nPntElem> &points, 
        Elastic *elas);
    ~SolidElement();
    
    // compute stiffness term
    void computeStiff() const;
    
    // measure cost 
    double measure(int count) const;
    
    // test stiffness 
    void test() const;
    
    // compute Real displacement, used by receiver
    void computeGroundMotion(Real phi, const RMatPP &weights, RRow3 &u_spz) const; 
    void computeStrain(Real phi, const RMatPP &weights, RRow6 &strain) const; 
    void computeCurl(Real phi, const RMatPP &weights, RRow3 &curl) const; 
    void forceTIso();
    
    // side-wise
    void feedDispOnSide(int side, CMatXX_RM &buffer, int row) const; 
    
    // verbose
    std::string verbose() const;
    
    // reset
    void resetZero();
    
private:
    
    // displ ==> stiff
    void displToStiff() const;
    
    // material
    Elastic *mElastic;
    CrdTransTIsoSolid *mCrdTransTIso;
    // flags
    bool mInTIso;
    bool mElem3D;
    
//-------------------------- static --------------------------//    
public:
    // initialize static workspace
    static void initWorkspace(int maxMaxNu);
    
private:
    static SolidResponse sResponse;
};


