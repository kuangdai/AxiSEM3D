// SFCoupling3D.cpp
// created by Kuangdai on 5-Apr-2016 
// 3D solid-fluid boundary condition

#include "SFCoupling3D.h"
#include "SolverFFTW_1.h"
#include "SolverFFTW_3.h"

void SFCoupling3D::coupleFluidToSolid(const CColX &fluidStiff, CMatX3 &solidStiff) const {
    // constants
    int Nr = mNormal_assembled_invMassFluid.rows();

    // copy
    CColX &stiffFluidC = SolverFFTW_1::getC2R_CMat(Nr);
    stiffFluidC = fluidStiff;
    
    // FFT forward    
    SolverFFTW_1::computeC2R(Nr);
    RColX &stiffFluidR = SolverFFTW_1::getC2R_RMat(Nr);
    RMatX3 &stiffSolidR = SolverFFTW_3::getR2C_RMat(Nr);
    
    // fluid => solid
    stiffSolidR.col(0) = mNormal_assembled_invMassFluid.col(0).schur(stiffFluidR);
    stiffSolidR.col(1) = mNormal_assembled_invMassFluid.col(1).schur(stiffFluidR);
    stiffSolidR.col(2) = mNormal_assembled_invMassFluid.col(2).schur(stiffFluidR);
    
    // FFT backward
    SolverFFTW_3::computeR2C(Nr);
    solidStiff -= SolverFFTW_3::getR2C_CMat(Nr);
}

void SFCoupling3D::coupleSolidToFluid(const CMatX3 &solidDispl, CColX &fluidStiff) const {
    // constants
    int Nr = mNormal_unassembled.rows();

    // copy
    CMatX3 &displSolidC = SolverFFTW_3::getC2R_CMat(Nr);
    displSolidC = solidDispl;
    
    // FFT forward    
    SolverFFTW_3::computeC2R(Nr);
    RMatX3 &displSolidR = SolverFFTW_3::getC2R_RMat(Nr);
    RColX &stiffFluidR = SolverFFTW_1::getR2C_RMat(Nr);
    
    // solid => fluid
    stiffFluidR = mNormal_unassembled.col(0).schur(displSolidR.col(0))
                + mNormal_unassembled.col(1).schur(displSolidR.col(1))
                + mNormal_unassembled.col(2).schur(displSolidR.col(2));
    
    // FFT backward
    SolverFFTW_1::computeR2C(Nr);
    fluidStiff += SolverFFTW_1::getR2C_CMat(Nr);
}


void SFCoupling3D::checkCompatibility(int nr) const {
    int myNr0 = mNormal_unassembled.rows();
    int myNr1 = mNormal_assembled_invMassFluid.rows();
    if (myNr0 != nr || myNr1 != nr) {
        throw std::runtime_error("SFCoupling3D::checkCompatibility || Incompatible size.");
    }
}

