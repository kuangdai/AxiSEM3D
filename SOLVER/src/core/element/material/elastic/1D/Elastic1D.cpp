// Elastic1D.h
// created by Kuangdai on 29-Apr-2016 
// base class of 1D elasticity

#include "Elastic1D.h"
#include "Attenuation1D.h"

Elastic1D::Elastic1D(Attenuation1D *att):
mAttenuation(att) {
    // nothing
}

Elastic1D::~Elastic1D() {
    if (mAttenuation) {
        delete mAttenuation;
    }
}

void Elastic1D::checkCompatibility(int Nr) const {
    if (mAttenuation) {
        mAttenuation->checkCompatibility(Nr);
    }
}

void Elastic1D::resetZero() {
    if (mAttenuation) {
        mAttenuation->resetZero();
    }
}

