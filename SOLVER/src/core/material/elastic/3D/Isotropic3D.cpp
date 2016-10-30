// Isotropic3D.cpp
// created by Kuangdai on 22-Apr-2016 
// isotropic 3D material

#include "Isotropic3D.h"
#include "SolverFFTW_N6.h"
#include "Attenuation3D.h"

Isotropic3D::Isotropic3D(const RMatXN &lambda, const RMatXN &mu, Attenuation3D *att):
Elastic3D(att), mLambda(lambda), mMu(mu), mMu2(two * mu) {
    // nothing
}

void Isotropic3D::strainToStress(const vec_ar9_CMatPP &strain, vec_ar9_CMatPP &stress, int dummy) const {
    // constants
    int Nr = mLambda.rows();
    int Nu = Nr / 2;
    
    // copy
    CMatXN6 &strainC = SolverFFTW_N6::getC2R_CMat(Nr);
    Elastic3D::flattenVectorVoigt(strain, strainC, Nu);
    
    // FFT forward
    SolverFFTW_N6::computeC2R(Nr);
    RMatXN6 &strainR = SolverFFTW_N6::getC2R_RMat(Nr);
    RMatXN6 &stressR = SolverFFTW_N6::getR2C_RMat(Nr);
    
    // strain => stress
    // to avoid dynamic allocation, use stressR.block(0, 3 * nPE, Nr, nPE) to store Sii
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
        
    // FFT backward
    SolverFFTW_N6::computeR2C(Nr);
    CMatXN6 &stressC = SolverFFTW_N6::getR2C_CMat(Nr);
    
    // copy
    Elastic3D::stackupVectorVoigt(stressC, stress, Nu);
}

void Isotropic3D::checkCompatibility(int Nr, bool isVoigt) const {
    Elastic3D::checkCompatibility(Nr, isVoigt);
    int myNr = mLambda.rows();
    if (Nr != myNr) 
        throw std::runtime_error("Isotropic3D::checkCompatibility || Incompatible size.");
    if (!isVoigt) throw std::runtime_error("Isotropic3D::checkCompatibility || "
        "Incompatible gradient operator.");
}


