// FluidPoint.h
// created by Kuangdai on 4-Apr-2016 
// fluid gll points 

#pragma once

class Mass;
#include "Point.h"

class FluidPoint: public Point {
    friend class SolidFluidPoint;

public:    
    FluidPoint(int nr, bool axial, const RDCol2 &crds, Mass *mass);
    ~FluidPoint();

    // update in time domain by Newmark
    void updateNewmark(Real dt);
    
    // check stability
    bool stable() const {return mDispl.allFinite();};
    
    // reset to zero 
    void resetZero(); 
    
    // randomize disp and stiff
    void randomDispl(Real factor = one, int seed = -1);
    void randomStiff(Real factor = one, int seed = -1);
    
    // verbose
    std::string verbose() const;
    
    // measure cost 
    double measure(int count);
    
    // test mass
    void test();
    
    // communication size
    int sizeComm() const {return mStiff.rows();};
    
    // communication
    void feedBuffer(CColX &buffer, int &row);
    void extractBuffer(CColX &buffer, int &row);
    
    ///////////// fluid-only /////////////
    // scatter displ to element
    void scatterDisplToElement(vec_CMatPP &displ, int ipol, int jpol, int maxNu) const;
    
    // gather stiff from element
    void gatherStiffFromElement(const vec_CMatPP &stiff, int ipol, int jpol);
    
    // wisdom
    void learnWisdom(Real cutoff);
    int getNuWisdom() const {return mNuWisdom;};
	
	// get displacement
	const CColX &getDispFourierFluid() const {return mDispl;};
    
private:
    
    // mask 
    void maskField(CColX &field);

    // fields
    CColX mDispl;
    CColX mVeloc;
    CColX mAccel;
    CColX mStiff;
    
    // mass
    Mass *mMass;
    
    // wisdom
    Real mMaxDisplWisdom = -1.;
    int mNuWisdom = 0;
};
