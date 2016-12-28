// SeismometerENZ.cpp
// created by Kuangdai on 8-Apr-2016 
// compute ground motion at a point in east, north, vertical

#include "SeismometerENZ.h"

SeismometerENZ::SeismometerENZ(Real phi, const RMatPP &weights, Element *element,
    Real theta, Real baz):
Seismometer(phi, weights, element), mTheta(theta), mBAz(baz) {
    // nothing
}

void SeismometerENZ::getGroundMotion(RRow3 &gm) const {
    Seismometer::getGroundMotion(gm);
    // first to r, theta, phi
    Real cost = cos(mTheta);
    Real sint = sin(mTheta);
    Real ur = gm(0) * sint + gm(2) * cost; 
    Real ut = gm(0) * cost - gm(2) * sint;
    Real up = gm(1);
    // then to east, north, vertical
    Real cosbaz = cos(mBAz);
    Real sinbaz = sin(mBAz);
    gm(0) = - ut * sinbaz + up * cosbaz;
    gm(1) = - ut * cosbaz - up * sinbaz;
    gm(2) = ur; 
}