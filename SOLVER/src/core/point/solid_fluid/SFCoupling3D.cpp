// SFCoupling3D.cpp
// created by Kuangdai on 5-Apr-2016 
// 3D solid-fluid boundary condition

#include "SFCoupling3D.h"
#include "SolverFFTW_1.h"
#include "SolverFFTW_3.h"

void SFCoupling3D::coupleFluidToSolid(const CColX &fluidStiff, CMatX3 &solidStiff) const {
    // constants
    int Nr = mNormal_assembled_invMassFluid.rows();
    int Nc = Nr / 2 + 1;
    
    // copy
    CColX &stiffFluidC = SolverFFTW_1::getC2R_CMat();
    stiffFluidC.topRows(Nc) = fluidStiff;
    
    // FFT forward    
    SolverFFTW_1::computeC2R(Nr);
    RColX &stiffFluidR = SolverFFTW_1::getC2R_RMat();
    RMatX3 &stiffSolidR = SolverFFTW_3::getR2C_RMat();
    
    // fluid => solid
    stiffSolidR.block(0, 0, Nr, 1) = mNormal_assembled_invMassFluid.col(0).schur(stiffFluidR.topRows(Nr));
    stiffSolidR.block(0, 1, Nr, 1) = mNormal_assembled_invMassFluid.col(1).schur(stiffFluidR.topRows(Nr));
    stiffSolidR.block(0, 2, Nr, 1) = mNormal_assembled_invMassFluid.col(2).schur(stiffFluidR.topRows(Nr));
    
    // FFT backward
    SolverFFTW_3::computeR2C(Nr);
    solidStiff -= SolverFFTW_3::getR2C_CMat().topRows(Nc);
}

void SFCoupling3D::coupleSolidToFluid(const CMatX3 &solidDispl, CColX &fluidStiff) const {
    // constants
    int Nr = mNormal_unassembled.rows();
    int Nc = Nr / 2 + 1;

    // copy
    CMatX3 &displSolidC = SolverFFTW_3::getC2R_CMat();
    displSolidC.topRows(Nc) = solidDispl;
    
    // FFT forward    
    SolverFFTW_3::computeC2R(Nr);
    RMatX3 &displSolidR = SolverFFTW_3::getC2R_RMat();
    RColX &stiffFluidR = SolverFFTW_1::getR2C_RMat();
    
    // solid => fluid
    stiffFluidR.topRows(Nr) = mNormal_unassembled.col(0).schur(displSolidR.block(0, 0, Nr, 1))
                            + mNormal_unassembled.col(1).schur(displSolidR.block(0, 1, Nr, 1))
                            + mNormal_unassembled.col(2).schur(displSolidR.block(0, 2, Nr, 1));
    
    // FFT backward
    SolverFFTW_1::computeR2C(Nr);
    fluidStiff += SolverFFTW_1::getR2C_CMat().topRows(Nc);
}


void SFCoupling3D::checkCompatibility(int nr) const {
    int myNr0 = mNormal_unassembled.rows();
    int myNr1 = mNormal_assembled_invMassFluid.rows();
    if (myNr0 != nr || myNr1 != nr) {
        throw std::runtime_error("SFCoupling3D::checkCompatibility || Incompatible size.");
    }
}

