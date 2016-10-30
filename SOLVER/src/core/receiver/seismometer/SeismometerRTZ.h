// SeismometerRTZ.h
// created by Kuangdai on 8-Apr-2016 
// compute ground motion at a point in radial, transverse, vertical

#pragma once

#include "Seismometer.h"

class SeismometerRTZ: public Seismometer {
public:    
    SeismometerRTZ(Real phi, const RMatPP &weights, Element *element, Real theta);

    void getGroundMotion(RRow3 &gm) const; 
        
protected:
    // distance
    Real mTheta;
};

