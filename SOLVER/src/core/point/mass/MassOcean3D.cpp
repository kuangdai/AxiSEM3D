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

    // copy 
    CMatX3 &stiffC = SolverFFTW_3::getC2R_CMat(Nr);
    stiffC = stiff;
    
    // FFT forward    
    SolverFFTW_3::computeC2R(Nr);
    RMatX3 &stiffR = SolverFFTW_3::getC2R_RMat(Nr);
    RMatX3 &accelR = SolverFFTW_3::getR2C_RMat(Nr);
    
    // stiff => accel
    accelR.col(0) = stiffR.col(0).schur(mInvMass);
    accelR.col(1) = stiffR.col(1).schur(mInvMass);
    accelR.col(2) = stiffR.col(2).schur(mInvMass);
    
    RColX &fn = SolverFFTW_1::getR2C_RMat(Nr);
    fn = stiffR.col(0).schur(mNormal_scal.col(0)) 
       + stiffR.col(1).schur(mNormal_scal.col(1))
       + stiffR.col(2).schur(mNormal_scal.col(2));
    
    accelR.col(0) -= fn.schur(mNormal_scal.col(0));
    accelR.col(1) -= fn.schur(mNormal_scal.col(1));
    accelR.col(2) -= fn.schur(mNormal_scal.col(2));
    
    // FFT backward
    SolverFFTW_3::computeR2C(Nr);
    stiff = SolverFFTW_3::getR2C_CMat(Nr);
}

void MassOcean3D::computeAccel(CColX &stiff) const {
    // constants
    int Nr = mInvMass.rows();

    // copy
    CColX &stiffC = SolverFFTW_1::getC2R_CMat(Nr);
    stiffC = stiff;
    
    // FFT forward    
    SolverFFTW_1::computeC2R(Nr);
    RColX &stiffR = SolverFFTW_1::getC2R_RMat(Nr);
    RColX &accelR = SolverFFTW_1::getR2C_RMat(Nr);
    
    // stiff => accel
    accelR = stiffR.schur(mInvMass);
    
    // FFT backward
    SolverFFTW_1::computeR2C(Nr);
    stiff = SolverFFTW_1::getR2C_CMat(Nr);
}

void MassOcean3D::checkCompatibility(int nr) const {
    if (mInvMass.rows() != nr) {
        throw std::runtime_error("MassOcean3D::checkCompatibility || Incompatible size.");
    }
}


