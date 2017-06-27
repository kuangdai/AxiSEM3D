// Volumetric3D_EMC.h
// created by Kuangdai on 16-May-2017 
// genetral Volumetric3D model with IRIS-EMC format

#pragma once

#include "Volumetric3D.h"
#include "eigenp.h"

class Volumetric3D_EMC: public Volumetric3D {
public:

    void initialize();
    void initialize(const std::vector<std::string> &params);
    bool get3dProperties(double r, double theta, double phi, double rElemCenter,
        std::vector<MaterialProperty> &properties, 
        std::vector<MaterialRefType> &refTypes,
        std::vector<double> &values) const;
    std::string verbose() const;
    
private:
    
    // file
    std::string mFileName;
    std::string mVarName;
    
    // property
    MaterialProperty mMaterialProp;
    MaterialRefType mReferenceType;
    
    // factor
    double mFactor = 1.0;
    
    // use geocentric or geographic
    bool mGeographic = false;
    
    // one-file-per-depth format
    bool mOneFilePerDepth = false;
    
    // data
    std::vector<RDMatXX> mGridData;
    RDColX mGridDep;
    RDColX mGridLat;
    RDColX mGridLon;
};

