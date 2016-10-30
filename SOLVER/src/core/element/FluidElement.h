// FluidElement.h
// created by Kuangdai on 29-Mar-2016 
// base class of fluid elements

#pragma once

#include "Element.h"

class Acoustic;

class FluidElement : public Element {
public:
    
    FluidElement(Gradient *grad, const std::array<Point *, nPntElem> &points, Acoustic *acous);
    ~FluidElement();
    
    // compute stiffness term
    void computeStiff() const;
    
    // measure cost
    double measure(int count, bool user) const;
    
    // test stiffness 
    void test() const;
    
    // compute Real displacement, used by receiver
    void computeGroundMotion(Real phi, const RMatPP &weights, RRow3 &u_spz) const; 
    
    // verbose
    std::string verbose() const;
    
private:
    
    // displ ==> stiff
    void displToStiff(const vec_CMatPP &displ, vec_CMatPP &stiff) const;
    
    // material
    Acoustic *mAcoustic;
    
//-------------------------- static --------------------------//    
public:
    // initialize static workspace
    static void initWorkspace(int maxMaxNu);
    
private:
    // static workspaces
    static vec_CMatPP sDispl;
    static vec_CMatPP sStiff;
    static vec_ar3_CMatPP sStrain;
    static vec_ar3_CMatPP sStress;
};
