// Point.h
// created by Kuangdai on 3-Apr-2016 
// base class of gll points

#pragma once

#include "global.h"
#include "eigenc.h"
#include "eigenp.h"

class Point {
public:    
    Point(int nr, bool axial, const RDCol2 &crds);
    virtual ~Point() {};
    
    // update in time domain by Newmark
    virtual void updateNewmark(double dt) = 0;
    
    // check stability
    virtual bool stable() const = 0;
    
    // reset to zero 
    virtual void resetZero() = 0; 
    
    // randomize disp and stiff
    virtual void randomDispl(Real factor = one, int seed = -1, int max_order = -1) = 0;
    virtual void randomStiff(Real factor = one, int seed = -1, int max_order = -1) = 0;
    
    // verbose
    virtual std::string verbose() const = 0;
    
    // measure cost
    virtual double measure(int count) = 0;
    
    // test mass
    virtual void test() = 0;
    
    // communication size
    virtual int sizeComm() const = 0;
    
    // communication
    virtual void feedBuffer(CColX &buffer, int &row) = 0;
    virtual void extractBuffer(CColX &buffer, int &row) = 0;
    
    // scatter displ to element
    virtual void scatterDisplToElement(vec_ar3_CMatPP &displ, int ipol, int jpol, int maxNu) const;
    virtual void scatterDisplToElement(vec_CMatPP &displ, int ipol, int jpol, int maxNu) const;
    
    // gather stiff from element
    virtual void gatherStiffFromElement(const vec_ar3_CMatPP &stiff, int ipol, int jpol);
    virtual void gatherStiffFromElement(const vec_CMatPP &stiff, int ipol, int jpol);
    
    // add to stiff, used by source
    virtual void addToStiff(const CMatX3 &source);
    
    // n_u, n_r
    int getNr() const {return mNr;};
    int getNu() const {return mNu;};
    
    // axis
    bool axial() const {return mAxial;};
    
    // location
    const RDCol2 &getCoords() const {return mCoords;};
    
    // wisdom 
    virtual void learnWisdom(Real cutoff) = 0;
    virtual int getNuWisdom() const = 0;
    
    // signature for cost measurement
    std::string costSignature() const;
    
    // get displacement
    virtual const CMatX3 &getDispFourierSolid() const;
    virtual const CColX &getDispFourierFluid() const;

protected:
    
    // number of slices on ring
    int mNr;   
    
    // Fourier expansion order of solution
    // mNu = mNr / 2 
    int mNu;
    
    // axis
    bool mAxial;
    
    // location (s,z)
    // used only for movie and error report
    RDCol2 mCoords;
    
public:
    // domain tag, mainly for debug
    void setDomainTag(int tag) {mDomainTag = tag;};
    int getDomainTag() const {return mDomainTag;};
    
private:
    int mDomainTag = -1;
    
};