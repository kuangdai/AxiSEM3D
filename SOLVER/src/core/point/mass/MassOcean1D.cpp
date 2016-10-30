// MassOcean1D.cpp
// created by Kuangdai on 3-Apr-2016 
// 1D mass with ocean

#include "MassOcean1D.h"

MassOcean1D::MassOcean1D(double mass, double massOcean, double theta) {
    mInvMassZ = (Real)(1. / (mass + massOcean));
    mInvMassR = (Real)(1. / mass);
    mSint = (Real)sin(theta);
    mCost = (Real)cos(theta);
}

void MassOcean1D::computeAccel(CMatX3 &stiff) const {
    CColX stiff_Z = stiff.col(0) * mSint + stiff.col(2) * mCost;
    CColX stiff_R = stiff.col(0) * mCost - stiff.col(2) * mSint;
    stiff_Z *= mInvMassZ;
    stiff_R *= mInvMassR;
    stiff.col(0) = stiff_Z * mSint + stiff_R * mCost;
    stiff.col(2) = stiff_Z * mCost - stiff_R * mSint;
    stiff.col(1) *= mInvMassR;
}

void MassOcean1D::computeAccel(CColX &stiff) const {
    stiff *= mInvMassR; 
}