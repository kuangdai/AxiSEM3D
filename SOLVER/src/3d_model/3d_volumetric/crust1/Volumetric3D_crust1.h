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
    void initialize(const std::vector<std::string> &params);
    
    bool get3dProperties(double r, double theta, double phi, double rElemCenter,
        double &vpv, double &vph, double &vsv, double &vsh, double &rho) const;
    
    ReferenceTypes getReferenceType() const {return ReferenceTypes::Absolute;};
    
    std::string verbose() const;
    
    // set outer radius
    void setROuter(double router) {mRSurf = router;};
    
    // obtain more mesh information from ExodusModel
    void setupExodusModel(const ExodusModel *exModel);
    
private:
    
    size_t columnSurf() const {
        size_t colSurf = 5; // no ice, no sediment
        if (mIncludeIce) {
            colSurf = 1; // ice
        } else if (mIncludeSediment) {
            colSurf = 2; // sediment
        }
        return colSurf;
    };
    
    static void interpThetaPhi(double theta, double phi, size_t np, 
        std::vector<size_t> &ilat, std::vector<size_t> &ilon, 
        std::vector<double> &wlat, std::vector<double> &wlon);
    
private:
    // model constants
    static const size_t sNLayer;
    static const size_t sNLat;
    static const size_t sNLon;
    
    // radii of reference sphere
    double mRSurf = 6371000.0;
    double mRMoho = 6346600.0;
    
    // options
    // include ice or not
    bool mIncludeIce = false;
    // include sediment or not
    bool mIncludeSediment = true;
    // number of interpolation points 
    size_t mNPointInterp = 2;
    // use geocentric or geographic
    bool mGeographic = false;
    
    // element boundaries in mesh
    std::vector<double> mElementBoundaries;
    
    // original data mapped onto reference sphere
    RDMatXX mRl;
    RDMatXX mVp;
    RDMatXX mVs;
    RDMatXX mRh;
    
    // thickness weighted data mapped onto reference sphere
    RDColX mRlGLL;
    RDMatXX mVpGLL;
    RDMatXX mVsGLL;
    RDMatXX mRhGLL;
};
