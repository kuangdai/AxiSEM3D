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
    static const size_t sNLayer;
    static const size_t sNLat;
    static const size_t sNLon;
    
    // smoothening 
    size_t mGaussianOrder = 2;
    double mGaussianDev = .5;
    
    // interpolation
    size_t mNPointInterp = 2;
    // use geocentric or geographic
    bool mGeographic = true;
    
    // treat ice as water load
    bool mIncludeIceAsWater = false;
    
    // flag to benchmark with specfem
    bool mBenchmarkSPECFEM = false; 
    
    // depth at grid points 
    RDMatXX mDepth;
};

