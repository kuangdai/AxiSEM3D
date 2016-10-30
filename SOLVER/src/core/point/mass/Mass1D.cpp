// Mass1D.cpp
// created by Kuangdai on 3-Apr-2016 
// 1D mass

#include "Mass1D.h"

Mass1D::Mass1D(Real invMass):
mInvMass(invMass) {
    // nothing
}

void Mass1D::computeAccel(CMatX3 &stiff) const {
    stiff *= mInvMass; 
}

void Mass1D::computeAccel(CColX &stiff) const {
    stiff *= mInvMass; 
}