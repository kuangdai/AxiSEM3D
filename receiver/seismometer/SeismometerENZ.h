// SeismometerENZ.h
// created by Kuangdai on 8-Apr-2016 
// compute ground motion at a point in east, north, vertical

#pragma once

#include "Seismometer.h"

class SeismometerENZ: public Seismometer {
public:    
    SeismometerENZ(Real phi, const RMatPP &weights, Element *element,
        Real theta, Real baz);
    void getGroundMotion(RRow3 &gm) const; 
        
protected:
    // distance
    Real mTheta;
    // back azimuth
    Real mBAz; 
    
};

