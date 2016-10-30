// SolidElement.h
// created by Kuangdai on 29-Mar-2016 
// base class of solid element

#pragma once

#include "Element.h"

class Elastic;

class SolidElement : public Element {
public:
    
    SolidElement(Gradient *grad, const std::array<Point *, nPntElem> &points, Elastic *elas);
    ~SolidElement();
    
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
    void displToStiff(const vec_ar3_CMatPP &displ, vec_ar3_CMatPP &stiff) const;
    
    // material
    Elastic *mElastic;
    
//-------------------------- static --------------------------//    
public:
    // initialize static workspace
    static void initWorkspace(int maxMaxNu);
    
private:
    // static workspaces
    static vec_ar3_CMatPP sDispl;
    static vec_ar3_CMatPP sStiff;
    static vec_ar9_CMatPP sStrain;
    static vec_ar9_CMatPP sStress;
};
