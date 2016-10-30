// Ellipticity.h
// created by Kuangdai on 4-Jun-2016 
// ellipticity

#pragma once

#include "Geometric3D.h"

class Ellipticity: public Geometric3D {
public:

    double getDeltaR(double r, double theta, double phi, double rElemCenter) const;
    
    // bool getNablaDeltaR(double r, double theta, double phi, double rElemCenter,
    //     double &deltaR_r, double &deltaR_theta, double &deltaR_phi) const;
    
    std::string verbose() const;
};

