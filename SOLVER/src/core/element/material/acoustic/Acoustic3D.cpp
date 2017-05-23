// Acoustic3D.cpp
// created by Kuangdai on 23-Apr-2016 
// 3D acoustic

#include "Acoustic3D.h"
#include "FluidElement.h"
#include "SolverFFTW_N3.h"

void Acoustic3D::strainToStress(FluidResponse &response) const {
    int Nr = response.mNr;
    const RMatXN3 &strain = SolverFFTW_N3::getC2R_RMat(Nr);
    RMatXN3 &stress = SolverFFTW_N3::getR2C_RMat(Nr);
    stress.block(0, 0 * nPE, Nr, nPE) = mKFlat.schur(strain.block(0, 0 * nPE, Nr, nPE));
    stress.block(0, 1 * nPE, Nr, nPE) = mKFlat.schur(strain.block(0, 1 * nPE, Nr, nPE)); 
    stress.block(0, 2 * nPE, Nr, nPE) = mKFlat.schur(strain.block(0, 2 * nPE, Nr, nPE));
}

void Acoustic3D::checkCompatibility(int Nr) const {
    int myNr = mKFlat.rows();
    if (Nr != myNr) {
        throw std::runtime_error("Acoustic3D::checkCompatibility || Incompatible size.");
    }
}
