// Isotropic3D.cpp
// created by Kuangdai on 22-Apr-2016 
// isotropic 3D material

#include "Isotropic3D.h"
#include "Attenuation3D.h"
#include "SolidElement.h"
#include "SolverFFTW_N6.h"

void Isotropic3D::strainToStress(SolidResponse &response) const {
    int Nr = response.mNr;
    // to avoid dynamic allocation, use stressR.block(0, 3 * nPE, Nr, nPE) to store Sii
    const RMatXN6 &strainR = SolverFFTW_N6::getC2R_RMat(Nr);
    RMatXN6 &stressR = SolverFFTW_N6::getR2C_RMat(Nr);
    stressR.block(0, 3 * nPE, Nr, nPE) = mLambda.schur(strainR.block(0, 0 * nPE, Nr, nPE) 
                                                     + strainR.block(0, 1 * nPE, Nr, nPE) 
                                                     + strainR.block(0, 2 * nPE, Nr, nPE));
    stressR.block(0, 0 * nPE, Nr, nPE) = stressR.block(0, 3 * nPE, Nr, nPE) + mMu2.schur(strainR.block(0, 0 * nPE, Nr, nPE));
    stressR.block(0, 1 * nPE, Nr, nPE) = stressR.block(0, 3 * nPE, Nr, nPE) + mMu2.schur(strainR.block(0, 1 * nPE, Nr, nPE));
    stressR.block(0, 2 * nPE, Nr, nPE) = stressR.block(0, 3 * nPE, Nr, nPE) + mMu2.schur(strainR.block(0, 2 * nPE, Nr, nPE));
    stressR.block(0, 3 * nPE, Nr, nPE) = mMu.schur(strainR.block(0, 3 * nPE, Nr, nPE));
    stressR.block(0, 4 * nPE, Nr, nPE) = mMu.schur(strainR.block(0, 4 * nPE, Nr, nPE));
    stressR.block(0, 5 * nPE, Nr, nPE) = mMu.schur(strainR.block(0, 5 * nPE, Nr, nPE));
    if (mAttenuation) {
        mAttenuation->applyToStress(stressR);
        mAttenuation->updateMemoryVariables(strainR);
    }
}

void Isotropic3D::checkCompatibility(int Nr) const {
    Elastic3D::checkCompatibility(Nr);
    int myNr = mLambda.rows();
    if (Nr != myNr) {
        throw std::runtime_error("Isotropic3D::checkCompatibility || Incompatible size.");
    } 
}

