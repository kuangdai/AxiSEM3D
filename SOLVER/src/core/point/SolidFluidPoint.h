// SolidFluidPoint.h
// created by Kuangdai on 5-Apr-2016 
// solid-fluid boundary condition

#pragma once

#include "Point.h"

class SolidPoint;
class FluidPoint;
class SFCoupling;

class SolidFluidPoint: public Point {
public:   
    SolidFluidPoint(SolidPoint *sp, FluidPoint *fp, SFCoupling *couple);
    ~SolidFluidPoint();
    
    void updateNewmark(double dt);
    
    // check stability
    bool stable() const;
    
    // reset to zero 
    void resetZero(); 
    
    // randomize disp and stiff
    void randomDispl(Real factor = one, int seed = -1, int max_order = -1);
    void randomStiff(Real factor = one, int seed = -1, int max_order = -1);
    
    // verbose
    std::string verbose() const;
    
    // measure cost
    double measure(int count);
    
    // test mass
    void test();
    
    // communication size
    int sizeComm() const;
    
    // communication
    void feedBuffer(CColX &buffer, int &row);
    void extractBuffer(CColX &buffer, int &row);
    
    // scatter displ to element
    void scatterDisplToElement(vec_ar3_CMatPP &displ, int ipol, int jpol, int maxNu) const;
    void scatterDisplToElement(vec_CMatPP &displ, int ipol, int jpol, int maxNu) const;
    
    // gather stiff from element
    void gatherStiffFromElement(const vec_ar3_CMatPP &stiff, int ipol, int jpol);
    void gatherStiffFromElement(const vec_CMatPP &stiff, int ipol, int jpol);
    
    // add to stiff, used by source
    void addToStiff(const CMatX3 &source);
    
    ///////////// solid-fluid-only /////////////   
    void coupleSolidFluid();
    
    // wisdom
    void learnWisdom(Real cutoff);
    int getNuWisdom() const;
    
    // get displacement
    const CMatX3 &getDispFourierSolid() const;
    const CColX &getDispFourierFluid() const;
    
private:
    double measureCoupling(int count);
    
    SolidPoint *mSolidPoint;
    FluidPoint *mFluidPoint;
    SFCoupling *mSFCoupling;
};
