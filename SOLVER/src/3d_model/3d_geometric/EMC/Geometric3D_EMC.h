// Geometric3D_EMC.h
// created by Kuangdai on 16-May-2017 
// topography on a boundary at any depth, with IRIS-EMC format

#pragma once

#include "Geometric3D.h"
#include "eigenp.h"

class Geometric3D_EMC: public Geometric3D {
public:

    void initialize();
    void initialize(const std::vector<std::string> &params);
    double getDeltaR(double r, double theta, double phi, double rElemCenter) const;
    std::string verbose() const;
    
private:
    
    // radius
    double mRLayer;
    double mRLower;
    double mRUpper;
    
    // file
    std::string mFileName;
    std::string mVarName;
    
    // factor
    double mFactor = 1.0;
    
    // use geocentric or geographic
    bool mGeographic = false;
    
    // data
    RDMatXX mGridData;
    RDColX mGridLat;
    RDColX mGridLon;
};

