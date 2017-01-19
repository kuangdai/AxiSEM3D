// Geometric3D_Internal.h
// created by Kuangdai on 19-Jan-2017 
// topography on a general internal boundary

#pragma once

#include "Geometric3D.h"
#include "eigenp.h"

class Geometric3D_Internal: public Geometric3D {
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
    
    // resolution
    int mN360;
    
    // file
    std::string mFileName;
    
    // use geocentric or geographic
    bool mGeographic = false;
    
    // factor
    double mFactor = 1.0;
    
    // smoothening 
    int mGaussianOrder = 0;
    double mGaussianDev = .5;

    // data
    RDMatXX mData;    
};

