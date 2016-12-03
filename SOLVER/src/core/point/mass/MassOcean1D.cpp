// MassOcean1D.cpp
// created by Kuangdai on 3-Apr-2016 
// 1D mass with ocean

#include "MassOcean1D.h"
#include "SolverFFTW_1.h"

MassOcean1D::MassOcean1D(double mass, double massOcean, double theta) {
    mInvMassZ = (Real)(1. / (mass + massOcean));
    mInvMassR = (Real)(1. / mass);
    mSint = (Real)sin(theta);
    mCost = (Real)cos(theta);
}

void MassOcean1D::computeAccel(CMatX3 &stiff) const {
    int nr_small = (stiff.rows() - 1) * 2;
    CColX &stiff_Z = SolverFFTW_1::getC2R_CMat(nr_small);
    CColX &stiff_R = SolverFFTW_1::getR2C_CMat(nr_small);
    stiff_Z = stiff.col(0) * mSint + stiff.col(2) * mCost;
    stiff_R = stiff.col(0) * mCost - stiff.col(2) * mSint;
    stiff_Z *= mInvMassZ;
    stiff_R *= mInvMassR;
    stiff.col(0) = stiff_Z * mSint + stiff_R * mCost;
    stiff.col(2) = stiff_Z * mCost - stiff_R * mSint;
    stiff.col(1) *= mInvMassR;
}

void MassOcean1D::computeAccel(CColX &stiff) const {
    stiff *= mInvMassR; 
}