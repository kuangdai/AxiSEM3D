// FluidElement.h
// created by Kuangdai on 29-Mar-2016 
// base class of fluid elements

#pragma once

#include "Element.h"

class Acoustic;
class CrdTransTIsoFluid;

// static workspaces
struct FluidResponse {
    // disp
    vec_CMatPP mDispl;
    // strain
    vec_ar3_CMatPP mStrain;
    // stress
    vec_ar3_CMatPP mStress;
    // stiff
    vec_CMatPP mStiff;
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


class FluidElement : public Element {
public:
    
    FluidElement(Gradient *grad, PRT *prt, const std::array<Point *, nPntElem> &points, 
        Acoustic *acous);
    ~FluidElement();
    
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
    
private:
    
    // displ ==> stiff
    void displToStiff() const;
    
    // material
    Acoustic *mAcoustic;
    CrdTransTIsoFluid *mCrdTransTIso;
    
    // flags
    bool mInTIso;
    bool mElem3D;
    
//-------------------------- static --------------------------//    
public:
    // initialize static workspace
    static void initWorkspace(int maxMaxNu);
    
private:
    // static workspaces
    static FluidResponse sResponse;
};
