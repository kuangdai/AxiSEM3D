// Volumetric3D_bubble.h
// created by Kuangdai on 16-May-2016 
// a bubble-shaped heterogeneity

#pragma once
#include "Volumetric3D.h"
#include "eigenp.h"

class Volumetric3D_bubble: public Volumetric3D {
public:
    
    void initialize(const std::vector<std::string> &params);
    
    bool get3dProperties(double r, double theta, double phi, double rElemCenter,
        std::vector<MaterialProperty> &properties, 
        std::vector<MaterialRefType> &refTypes,
        std::vector<double> &values) const;
    
    std::string verbose() const;
    
    void setSourceLocation(double srcLat, double srcLon, double srcDep) {
        mSrcLat = srcLat;
        mSrcLon = srcLon;
        mSrcDep = srcDep;
    }
    
    bool makeFluid3D() const {return mFluid;};
    
private:
    // property name
    MaterialProperty mMaterialProp;
    
    // value inside the bubble
    double mValueInside; 

    // reference type
    MaterialRefType mReferenceType;
    
    // radius of the bubble
    double mRadius;

    // center of the bubble
    double mDepth;
    double mLat, mLon;
    
    // source-centered
    bool mSourceCentered = false;
    double mSrcLat = 0.;
    double mSrcLon = 0.;
    double mSrcDep = 0.;
    
    // 3d fluid
    bool mFluid = false;
    
    // halfwidth at half maximum of Gaussian 
    // how the perturbation fades outside the bubble
    double mHWHM = -1.;
    
    // temp variables for performance
    RDCol3 mXyzBubble;
};
