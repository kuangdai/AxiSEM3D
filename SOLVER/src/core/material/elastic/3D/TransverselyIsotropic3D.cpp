// TransverselyIsotropic3D.cpp
// created by Kuangdai on 22-Apr-2016 
// transversely isotropic 3D material

#include "TransverselyIsotropic3D.h"
#include "SolverFFTW_N6.h"
#include "Attenuation3D.h"

TransverselyIsotropic3D::TransverselyIsotropic3D(const RDRowN &theta, 
    const RMatXN &A, const RMatXN &C, const RMatXN &F, 
    const RMatXN &L, const RMatXN &N, Attenuation3D *att): 
Elastic3D(att), mA(A), mC(C), mF(F), mL(L), mN(N), mN2(two * N) {
    int Nr = mA.rows();
    int n = Nr / 2 + 1;
    mSin1t = (theta.array().sin().matrix()).cast<Real>().colwise().replicate(n);
    mCos1t = (theta.array().cos().matrix()).cast<Real>().colwise().replicate(n);
    mSin2t = ((2. * theta).array().sin().matrix()).cast<Real>().colwise().replicate(n);
    mCos2t = ((2. * theta).array().cos().matrix()).cast<Real>().colwise().replicate(n);
}

void TransverselyIsotropic3D::strainToStress(const vec_ar9_CMatPP &strain, vec_ar9_CMatPP &stress, int dummy) const {
    // constants
    int Nr =  mA.rows();
    int Nu = Nr / 2;
    
    // copy forward
    CMatXN6 &strainCopy = SolverFFTW_N6::getR2C_CMat(Nr);
    Elastic3D::flattenVectorVoigt(strain, strainCopy, Nu);
    
    // rotate forward
    CMatXN6 &strainTIsoC = SolverFFTW_N6::getC2R_CMat(Nr);
    rotateStrainToTIso(strainCopy, strainTIsoC, Nu);
    
    // FFT forward
    SolverFFTW_N6::computeC2R(Nr);
    RMatXN6 &strainTIsoR = SolverFFTW_N6::getC2R_RMat(Nr);
    RMatXN6 &stressTIsoR = SolverFFTW_N6::getR2C_RMat(Nr);
    
    // strain => stress
    // to avoid dynamic allocation, use stressTIsoR.block(0, 3&4 * nPE, Nr, nPE) as temp memory
    stressTIsoR.block(0, 3 * nPE, Nr, nPE) = strainTIsoR.block(0, 0 * nPE, Nr, nPE) + strainTIsoR.block(0, 1 * nPE, Nr, nPE);
    stressTIsoR.block(0, 4 * nPE, Nr, nPE) = mA.schur(stressTIsoR.block(0, 3 * nPE, Nr, nPE)) + mF.schur(strainTIsoR.block(0, 2 * nPE, Nr, nPE));
    
    stressTIsoR.block(0, 0 * nPE, Nr, nPE) = stressTIsoR.block(0, 4 * nPE, Nr, nPE) - mN2.schur(strainTIsoR.block(0, 1 * nPE, Nr, nPE));
    stressTIsoR.block(0, 1 * nPE, Nr, nPE) = stressTIsoR.block(0, 4 * nPE, Nr, nPE) - mN2.schur(strainTIsoR.block(0, 0 * nPE, Nr, nPE));
    stressTIsoR.block(0, 2 * nPE, Nr, nPE) = mC.schur(strainTIsoR.block(0, 2 * nPE, Nr, nPE)) + mF.schur(stressTIsoR.block(0, 3 * nPE, Nr, nPE));
    stressTIsoR.block(0, 3 * nPE, Nr, nPE) = mL.schur(strainTIsoR.block(0, 3 * nPE, Nr, nPE));
    stressTIsoR.block(0, 4 * nPE, Nr, nPE) = mL.schur(strainTIsoR.block(0, 4 * nPE, Nr, nPE));
    stressTIsoR.block(0, 5 * nPE, Nr, nPE) = mN.schur(strainTIsoR.block(0, 5 * nPE, Nr, nPE));
    
    if (mAttenuation) {
        mAttenuation->applyToStress(stressTIsoR);
        mAttenuation->updateMemoryVariables(strainTIsoR);
    }
        
    // FFT backward
    SolverFFTW_N6::computeR2C(Nr);
    CMatXN6 &stressTIsoC = SolverFFTW_N6::getR2C_CMat(Nr);
    
    // rotate backward
    CMatXN6 &stressCopy = SolverFFTW_N6::getC2R_CMat(Nr);
    rotateStressToCyln(stressTIsoC, stressCopy, Nu);
    
    // copy backward
    Elastic3D::stackupVectorVoigt(stressCopy, stress, Nu);
}

void TransverselyIsotropic3D::checkCompatibility(int Nr, bool isVoigt) const {
    Elastic3D::checkCompatibility(Nr, isVoigt);
    int myNr = mA.rows();
    if (Nr != myNr) 
        throw std::runtime_error("TransverselyIsotropic3D::checkCompatibility || Incompatible size.");
    if (!isVoigt) throw std::runtime_error("TransverselyIsotropic3D::checkCompatibility || "
        "Incompatible gradient operator.");
}


void TransverselyIsotropic3D::rotateStrainToTIso(const CMatXN6 &strainCyln, CMatXN6 &strainTIso, int Nu) const {
    int n = Nu + 1;
    // temp
    strainTIso.block(0, nPE * 3, n, nPE) = strainCyln.block(0, nPE * 0, n, nPE) + strainCyln.block(0, nPE * 2, n, nPE); 
    strainTIso.block(0, nPE * 5, n, nPE) = strainCyln.block(0, nPE * 0, n, nPE) - strainCyln.block(0, nPE * 2, n, nPE); 
    // second order
    strainTIso.block(0, nPE * 0, n, nPE) = half * (strainTIso.block(0, nPE * 3, n, nPE) + mCos2t.schur(strainTIso.block(0, nPE * 5, n, nPE)) - mSin2t.schur(strainCyln.block(0, nPE * 4, n, nPE)));
    strainTIso.block(0, nPE * 2, n, nPE) = strainTIso.block(0, nPE * 3, n, nPE) - strainTIso.block(0, nPE * 0, n, nPE);
    strainTIso.block(0, nPE * 4, n, nPE) = mCos2t.schur(strainCyln.block(0, nPE * 4, n, nPE)) + mSin2t.schur(strainTIso.block(0, nPE * 5, n, nPE));
    // first order
    strainTIso.block(0, nPE * 3, n, nPE) = mCos1t.schur(strainCyln.block(0, nPE * 3, n, nPE)) + mSin1t.schur(strainCyln.block(0, nPE * 5, n, nPE));
    strainTIso.block(0, nPE * 5, n, nPE) = mCos1t.schur(strainCyln.block(0, nPE * 5, n, nPE)) - mSin1t.schur(strainCyln.block(0, nPE * 3, n, nPE));
    // zeroth order
    strainTIso.block(0, nPE * 1, n, nPE) = strainCyln.block(0, nPE * 1, n, nPE);
}
    
void TransverselyIsotropic3D::rotateStressToCyln(const CMatXN6 &stressTIso, CMatXN6 &stressCyln, int Nu) const {
    int n = Nu + 1;
    // temp
    stressCyln.block(0, nPE * 3, n, nPE) = stressTIso.block(0, nPE * 0, n, nPE) + stressTIso.block(0, nPE * 2, n, nPE); 
    stressCyln.block(0, nPE * 5, n, nPE) = (stressTIso.block(0, nPE * 0, n, nPE) - stressTIso.block(0, nPE * 2, n, nPE)) * half; 
    // second order
    stressCyln.block(0, nPE * 0, n, nPE) = half * stressCyln.block(0, nPE * 3, n, nPE) + mCos2t.schur(stressCyln.block(0, nPE * 5, n, nPE)) + mSin2t.schur(stressTIso.block(0, nPE * 4, n, nPE));
    stressCyln.block(0, nPE * 2, n, nPE) = stressCyln.block(0, nPE * 3, n, nPE) - stressCyln.block(0, nPE * 0, n, nPE);
    stressCyln.block(0, nPE * 4, n, nPE) = mCos2t.schur(stressTIso.block(0, nPE * 4, n, nPE)) - mSin2t.schur(stressCyln.block(0, nPE * 5, n, nPE));
    // first order
    stressCyln.block(0, nPE * 3, n, nPE) = mCos1t.schur(stressTIso.block(0, nPE * 3, n, nPE)) - mSin1t.schur(stressTIso.block(0, nPE * 5, n, nPE));
    stressCyln.block(0, nPE * 5, n, nPE) = mCos1t.schur(stressTIso.block(0, nPE * 5, n, nPE)) + mSin1t.schur(stressTIso.block(0, nPE * 3, n, nPE));
    // zeroth order
    stressCyln.block(0, nPE * 1, n, nPE) = stressTIso.block(0, nPE * 1, n, nPE);
}    

