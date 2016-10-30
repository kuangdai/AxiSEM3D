// Seismometer.cpp
// created by Kuangdai on 5-Apr-2016 
// compute ground motion at a point
// implement a new sub-class for other trace components 

#include "Seismometer.h"
#include "Element.h"

Seismometer::Seismometer(Real phi, const RMatPP &weights, Element *element):
mPhi(phi), mWeights(weights), mElement(element) {
    // nothing
}

void Seismometer::getGroundMotion(RRow3 &gm) const {
    mElement->computeGroundMotion(mPhi, mWeights, gm);
}
