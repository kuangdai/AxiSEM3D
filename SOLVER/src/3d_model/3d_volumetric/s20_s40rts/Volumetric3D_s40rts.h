// Volumetric3D_s40rts.h
// created by Kuangdai on 16-May-2016 
// mantle model s40rts
// http://www.earth.lsa.umich.edu/~jritsema/research.html

#pragma once
#include "Volumetric3D.h"

class Volumetric3D_s40rts: public Volumetric3D {
public:
    
    void initialize();
    void initialize(const std::vector<std::string> &params);
    
    void finalize();
    
    bool get3dProperties(double r, double theta, double phi, double rElemCenter,
        double &dvpv, double &dvph, double &dvsv, double &dvsh, double &drho) const;
    
    ReferenceTypes getReferenceType() const {return ReferenceTypes::Reference1D;};
    
    std::string verbose() const;
    
    // set outer radius
    void setROuter(double router) {mRSurf = router;};
    
private:
    double mRCMB = 3480000.0;
    double mRMoho = 6346600.0;
    double mRSurf = 6371000.0;
    double mScaleRho = .4;
};
