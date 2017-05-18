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
        std::vector<MaterialProperty> &properties, 
        std::vector<MaterialRefType> &refTypes,
        std::vector<double> &values) const;
    
    std::string verbose() const;
    
private:
    double mRCMB = 3480000.0;
    double mRMoho = 6346600.0;
    double mRSurf = 6371000.0;
    double mScaleRho = .4;
};
