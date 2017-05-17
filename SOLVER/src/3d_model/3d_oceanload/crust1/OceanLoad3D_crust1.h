// OceanLoad3D_crust1.h
// created by Kuangdai on 8-Oct-2016 
// crustal model CRUST 1.0 
// http://igppweb.ucsd.edu/~gabi/crust1.html

#pragma once

#include "OceanLoad3D.h"
#include "eigenp.h"

class OceanLoad3D_crust1: public OceanLoad3D {
public:

    void initialize();
    void initialize(const std::vector<std::string> &params);
    
    double getOceanDepth(double theta, double phi) const;
    
    std::string verbose() const;
    
private:
    
    // model constants
    static const int sNLayer;
    static const int sNLat;
    static const int sNLon;
    
    // smoothening 
    int mGaussianOrder = 2;
    double mGaussianDev = .5;
    
    // use geocentric or geographic
    bool mGeographic = true;
    
    // treat ice as water load
    bool mIncludeIceAsWater = false;
    
    // flag to benchmark with specfem
    bool mBenchmarkSPECFEM = false; 
    
    // depth at grid points 
    RDMatXX mDepth;
    RDColX mGridLat, mGridLon;
};

