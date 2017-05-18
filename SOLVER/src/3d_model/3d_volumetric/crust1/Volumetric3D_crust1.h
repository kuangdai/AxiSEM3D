// Volumetric3D_crust1.h
// created by Kuangdai on 16-May-2016 
// crustal model CRUST 1.0 
// http://igppweb.ucsd.edu/~gabi/crust1.html

#pragma once
#include "Volumetric3D.h"
#include "eigenp.h"

class Volumetric3D_crust1: public Volumetric3D {
    
public:
    void initialize();
    void initialize(const std::vector<std::string> &params);
    
    bool get3dProperties(double r, double theta, double phi, double rElemCenter,
        std::vector<MaterialProperty> &properties, 
        std::vector<MaterialRefType> &refTypes,
        std::vector<double> &values) const;
    
    std::string verbose() const;
    
    void setupExodusModel(const ExodusModel *exModel);
    
private:
    
    int columnSurf() const {
        int colSurf = 5; // no ice, no sediment
        if (mIncludeIce) {
            colSurf = 1; // ice
        } else if (mIncludeSediment) {
            colSurf = 2; // sediment
        }
        return colSurf;
    };
    
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
    // include ice or not
    bool mIncludeIce = false;
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
    
    // lat and lon grid
    RDColX mGridLat, mGridLon;
};
