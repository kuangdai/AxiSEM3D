// Seismometer.h
// created by Kuangdai on 5-Apr-2016 
// compute ground motion at a point
// implement a new sub-class for other trace components 

#pragma once

#include "eigenc.h"

class Element;

class Seismometer {
public:    
    Seismometer(Real phi, const RMatPP &weights, Element *element);
    virtual ~Seismometer() {};
    
    virtual void getGroundMotion(RRow3 &gm) const; 
    
protected:
    // azimuth
    Real mPhi;
    // interpolation weights
    RMatPP mWeights;
    // host Element
    Element *mElement;
};