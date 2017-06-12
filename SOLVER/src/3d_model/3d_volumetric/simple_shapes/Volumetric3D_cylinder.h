// Volumetric3D_cylinder.h
// created by Kuangdai on 16-May-2016 
// a cylinder-shaped heterogeneity

#pragma once
#include "Volumetric3D.h"
#include "eigenp.h"

class Volumetric3D_cylinder: public Volumetric3D {
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
    
    // value inside the cylinder
    double mValueInside; 

    // reference type
    MaterialRefType mReferenceType;
    
    // radius of the cylinder
    double mRadius;

    // anchor points of the cylinder
    double mD1, mD2;
    double mLat1, mLat2;
    double mLon1, mLon2;
    
    // source-centered
    bool mSourceCentered = false;
    double mSrcLat = 0.;
    double mSrcLon = 0.;
    double mSrcDep = 0.;
    
    // 3d fluid
    bool mFluid = false;
    
    // halfwidth at half maximum of Gaussian 
    // how the perturbation fades laterally and longitudinally outside the cylinder
    double mHWHM_lateral = -1.;
    double mHWHM_top_bot = -1.; 
    
    // temp variables for performance
    RDCol3 mXyzPoint1, mXyzPoint2;
    double mLength;
};
