// Geometric3D_crust1.h
// created by Kuangdai on 17-Jun-2016 
// crustal model CRUST 1.0 
// http://igppweb.ucsd.edu/~gabi/crust1.html


#pragma once

#include "Geometric3D.h"
#include "eigenp.h"

class Geometric3D_crust1: public Geometric3D {
public:

    void initialize();
    void initialize(const std::vector<double> &params);
    
    double getDeltaR(double r, double theta, double phi, double rElemCenter) const;
    
    // bool getNablaDeltaR(double r, double theta, double phi, double rElemCenter,
    //     double &deltaR_r, double &deltaR_theta, double &deltaR_phi) const;
    
    std::string verbose() const;
    
    // set outer radius
    void setROuter(double router) {mRSurf = router;};
    
private:
    
    // model constants
    static const int sNLayer;
    static const int sNLat;
    static const int sNLon;
    
    // radii of reference spheres
    double mRSurf = 6371000.0;
    double mRMoho = 6346600.0;
    double mRBase = 6151000.0; // the radius where deltaR is fixed at zero
    
    // include sediment or not
    bool mIncludeSediment = true;
    
    // strengthening factor
    double mSurfFactor = 1.;
    double mMohoFactor = 1.;
    
    // smoothening 
    int mGaussianOrder = 2;
    double mGaussianDev = .5;
    
    // interpolation
    int mNPointInterp = 2;
    // use geocentric or geographic
    bool mGeographic = true;
    
    // deltaR at surface and moho 
    RDMatXX mDeltaRSurf;
    RDMatXX mDeltaRMoho;
    
    // precomputed polar values
    // void computePolar();
    // RDCol2 mDrDxNorthSurf;
    // RDCol2 mDrDxSouthSurf;
    // RDCol2 mDrDxNorthMoho;
    // RDCol2 mDrDxSouthMoho;
};

