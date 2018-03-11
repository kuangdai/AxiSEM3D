// GLLPoint.h
// created by Kuangdai on 6-May-2016 
// general gll point 

#pragma once

#include "eigenp.h"

class Domain;

class GLLPoint {
public:
    // called from Mesh
    GLLPoint();
    
    // called from Quad 
    void setup(int nr, bool axial, bool surface, const RDCol2 &crds, double distTol);
    void addMassSolid(const RDColX &mass) {mMassSolid += mass;};
    void addMassFluid(const RDColX &mass) {mMassFluid += mass;};
    void addSFNormal(const RDMatX3 &normal) {
        mSFNormal += normal; 
        mSFNormal_assmble += normal;
    };
    
    void setOceanDepth(const RDColX &depth) {mOceanDepth = depth;};
    void addSurfNormal(const RDMatX3 &normal) {mSurfNormal += normal;};
    
    // release to domain
    int release(Domain &domain) const;
    
    // gets 
    int getNr() const {return mNr;};
    int getReferenceCount() const {return mReferenceCount;};
    const RDCol2 &getCoords() const {return mCoords;};
    
    // feed/extract buffer
    void feedBuffer(RDMatXX &buffer, int col);
    void extractBuffer(RDMatXX &buffer, int col);
    
private:
    // properties
    int mNr;
    bool mIsAxial;
    bool mOnSurface;
    RDCol2 mCoords;
    
    // mass
    RDColX mMassSolid;
    RDColX mMassFluid;
    
    // solid-fluid 
    RDMatX3 mSFNormal;
    RDMatX3 mSFNormal_assmble;
    
    // additional mass due to ocean load
    RDColX mOceanDepth;
    // surface normal 
    // NOTE: it should be impossible for one point to be located on both 
    // surface and s-f boundary. However, we don't want to mix the usage
    // of the two normal vectors. 
    RDMatX3 mSurfNormal;
    
    // total reference count
    int mReferenceCount;
    
};

