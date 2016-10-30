// SFCoupling1D.cpp
// created by Kuangdai on 5-Apr-2016 
// 1D solid-fluid boundary condition

#include "SFCoupling1D.h"
#include "SolidPoint.h"
#include "FluidPoint.h"

void SFCoupling1D::coupleFluidToSolid(const CColX &fluidStiff, CMatX3 &solidStiff) const {
    solidStiff.col(0) -= mNormalS_assembled_invMassFluid * fluidStiff;
    solidStiff.col(2) -= mNormalZ_assembled_invMassFluid * fluidStiff;
}

void SFCoupling1D::coupleSolidToFluid(const CMatX3 &solidDispl, CColX &fluidStiff) const {
    fluidStiff += mNormalS_unassembled * solidDispl.col(0) 
                + mNormalZ_unassembled * solidDispl.col(2);
}

