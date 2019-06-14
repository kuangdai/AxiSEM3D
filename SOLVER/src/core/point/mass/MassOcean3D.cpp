// MassOcean3D.cpp
// created by Kuangdai on 2-Jun-2016 
// 3D mass with ocean

#include "MassOcean3D.h"
#include "SolverFFTW_3.h"
#include "SolverFFTW_1.h"

MassOcean3D::MassOcean3D(const RDColX &mass, const RDColX &massOcean, const RDMatX3 &normal) {
    mInvMass = mass.array().pow(-1.).matrix().cast<Real>();
    RColX scal = (massOcean.array() / mass.schur(mass + massOcean).array()).sqrt().matrix().cast<Real>();
    mNormal_scal = normal.cast<Real>();
    mNormal_scal.col(0).array() *= scal.array();
    mNormal_scal.col(1).array() *= scal.array();
    mNormal_scal.col(2).array() *= scal.array();
}

void MassOcean3D::computeAccel(CMatX3 &stiff) const {
    // constants
    int Nr = mInvMass.rows();
    int Nc = Nr / 2 + 1;

    // copy 
    CMatX3 &stiffC = SolverFFTW_3::getC2R_CMat();
    stiffC.topRows(Nc) = stiff;
    
    // FFT forward    
    SolverFFTW_3::computeC2R(Nr);
    RMatX3 &stiffR = SolverFFTW_3::getC2R_RMat();
    RMatX3 &accelR = SolverFFTW_3::getR2C_RMat();
    
    // stiff => accel
    accelR.block(0, 0, Nr, 1) = stiffR.block(0, 0, Nr, 1).schur(mInvMass);
    accelR.block(0, 1, Nr, 1) = stiffR.block(0, 1, Nr, 1).schur(mInvMass);
    accelR.block(0, 2, Nr, 1) = stiffR.block(0, 2, Nr, 1).schur(mInvMass);
    
    RColX &fn = SolverFFTW_1::getR2C_RMat();
    fn.topRows(Nr) = stiffR.block(0, 0, Nr, 1).schur(mNormal_scal.col(0)) 
                   + stiffR.block(0, 1, Nr, 1).schur(mNormal_scal.col(1))
                   + stiffR.block(0, 2, Nr, 1).schur(mNormal_scal.col(2));
    
    accelR.block(0, 0, Nr, 1) -= fn.topRows(Nr).schur(mNormal_scal.col(0));
    accelR.block(0, 1, Nr, 1) -= fn.topRows(Nr).schur(mNormal_scal.col(1));
    accelR.block(0, 2, Nr, 1) -= fn.topRows(Nr).schur(mNormal_scal.col(2));
    
    // FFT backward
    SolverFFTW_3::computeR2C(Nr);
    stiff = SolverFFTW_3::getR2C_CMat().topRows(Nc);
}

void MassOcean3D::computeAccel(CColX &stiff) const {
    throw std::runtime_error("MassOcean3D::checkCompatibility || Fluid point with ocean load.");
}

void MassOcean3D::checkCompatibility(int nr) const {
    if (mInvMass.rows() != nr) {
        throw std::runtime_error("MassOcean3D::checkCompatibility || Incompatible size.");
    }
}


