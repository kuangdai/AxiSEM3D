// Mass3D.cpp
// created by Kuangdai on 2-Jun-2016 
// 3D mass

#include "Mass3D.h"
#include "SolverFFTW_3.h"
#include "SolverFFTW_1.h"

Mass3D::Mass3D(const RColX &invMass): mInvMass(invMass) {
    // nothing
}

void Mass3D::computeAccel(CMatX3 &stiff) const {
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
    
    // FFT backward
    SolverFFTW_3::computeR2C(Nr);
    stiff = SolverFFTW_3::getR2C_CMat(Nr);
}

void Mass3D::computeAccel(CColX &stiff) const {
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

void Mass3D::checkCompatibility(int nr) const {
    int myNr = mInvMass.rows();
    if (myNr != nr) {
        throw std::runtime_error("Mass3D::checkCompatibility || Incompatible size.");
    }
}


