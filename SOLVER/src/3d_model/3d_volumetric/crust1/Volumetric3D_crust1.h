// Volumetric3D_crust1.h
// created by Kuangdai on 16-May-2016 
// crustal model CRUST 1.0 
// http://igppweb.ucsd.edu/~gabi/crust1.html

#pragma once
#include "Volumetric3D.h"
#include "eigenp.h"

class Volumetric3D_crust1: public Volumetric3D {
    friend class Geometric3D_crust1;
    friend class OceanLoad3D_crust1;
    
public:
    
    void initialize();
    void initialize(const std::vector<double> &params);
    
    bool get3dProperties(double r, double theta, double phi, double rElemCenter,
        double &vpv, double &vph, double &vsv, double &vsh, double &rho) const;
    
    ReferenceTypes getReferenceType() const {return ReferenceTypes::Absolute;};
    
    std::string verbose() const;
    
    // set outer radius
    void setROuter(double router) {mRSurf = router;};
    
private:
    static void interpThetaPhi(double theta, double phi, int np, 
        std::vector<int> &ilat, std::vector<int> &ilon, 
        std::vector<double> &wlat, std::vector<double> &wlon);
    
private:
    // model constants
    static const int sNLayer;
    static const int sNLat;
    static const int sNLon;
    
    // radii of reference sphere
    double mRSurf = 6371000.0;
    double mRMoho = 6346600.0;
    
    // options
    // include sediment or not
    bool mIncludeSediment = true;
    double mMinimumSedimentThickness = 2000.0;
    // number of interpolation points 
    int mNPointInterp = 2;
    // use geocentric or geographic
    bool mGeographic = false;
    
    // data mapped onto reference sphere
    RDMatXX mRl;
    RDMatXX mVp;
    RDMatXX mVs;
    RDMatXX mRh;    
};
