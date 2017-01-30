// Volumetric3D_bubble.h
// created by Kuangdai on 16-May-2016 
// a bubble-shaped heterogeneity

#pragma once
#include "Volumetric3D.h"

class Volumetric3D_bubble: public Volumetric3D {
public:
    
    void initialize(const std::vector<std::string> &params);
    
    bool get3dProperties(double r, double theta, double phi, double rElemCenter,
        double &dvpv, double &dvph, double &dvsv, double &dvsh, double &drho) const;
    
    ReferenceTypes getReferenceType() const {return mReferenceType;};
    
    std::string verbose() const;
    
    void setSource(double srcLat, double srcLon, double srcDep) {
        mSrcLat = srcLat;
        mSrcLon = srcLon;
        mSrcDep = srcDep;
    }
    
private:
    // center of the bubble
    double mDepth;
    double mLat;
    double mLon;
    
    // Gaussian parameters
    double mRadius; // radius of the bubble 
    double mHWHM;   // halfwidth at half maximum; how the perturbation fades outside the bubble 
    double mMax;    // value inside the bubble
    
    // reference type
    ReferenceTypes mReferenceType;
    
    // source-centered
    bool mSourceCentered = false;
    double mSrcLat = 0.;
    double mSrcLon = 0.;
    double mSrcDep = 0.;
    
    // optional 
    bool mChangeVp = true;
    bool mChangeVs = true;
    bool mChangeRho = true;
};
