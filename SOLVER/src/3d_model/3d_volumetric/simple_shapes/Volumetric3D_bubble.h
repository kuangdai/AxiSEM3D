// Volumetric3D_bubble.h
// created by Kuangdai on 16-May-2016 
// a bubble-shaped heterogeneity

#pragma once
#include "Volumetric3D.h"

class Volumetric3D_bubble: public Volumetric3D {
public:
    
    void initialize(const std::vector<double> &params);
    
    bool get3dProperties(double r, double theta, double phi, double rElemCenter,
        double &dvpv, double &dvph, double &dvsv, double &dvsh, double &drho) const;
    
    ReferenceTypes getReferenceType() const {return mReferenceType;};
    
    std::string verbose() const;
    
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
    
    // optional 
    bool mChangeVp = true;
    bool mChangeVs = true;
    bool mChangeRho = true;
};
