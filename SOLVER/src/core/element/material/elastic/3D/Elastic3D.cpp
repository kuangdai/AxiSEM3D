// Elastic3D.h
// created by Kuangdai on 29-Apr-2016 
// base class of 3D elasticity 

#include "Elastic3D.h"
#include "Attenuation3D.h"

Elastic3D::Elastic3D(Attenuation3D *att):
mAttenuation(att) {
    // nothing
}

Elastic3D::~Elastic3D() {
    if (mAttenuation) {
        delete mAttenuation;
    }
}

void Elastic3D::checkCompatibility(int Nr) const {
    if (mAttenuation) {
        mAttenuation->checkCompatibility(Nr);
    }
}

void Elastic3D::resetZero() {
    if (mAttenuation) {
        mAttenuation->resetZero();
    }
}
