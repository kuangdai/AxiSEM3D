// SeismometerRTZ.cpp
// created by Kuangdai on 8-Apr-2016 
// compute ground motion at a point in radial, transverse, vertical

#include "SeismometerRTZ.h"

SeismometerRTZ::SeismometerRTZ(Real phi, const RMatPP &weights, Element *element, Real theta):
Seismometer(phi, weights, element), mTheta(theta) {
    // nothing
}

void SeismometerRTZ::getGroundMotion(RRow3 &gm) const {
    Seismometer::getGroundMotion(gm);
    Real cost = cos(mTheta);
    Real sint = sin(mTheta);
    Real u_s = gm(0);
    Real u_z = gm(2);
    gm(0) = u_s * cost - u_z * sint;
    // gm(1) = u_p, nothing to do
    gm(2) = u_s * sint + u_z * cost;  
}

