// TransverselyIsotropic3X.cpp
// created by Kuangdai on 6-Jun-2016 
// isotropic 3D material with topography

#include "TransverselyIsotropic3X.h"
#include "SolverFFTW_N9.h"
#include "SolverFFTW_N6.h"
#include "Attenuation3D.h"

TransverselyIsotropic3X::TransverselyIsotropic3X(const RDRowN &theta, 
    const RMatXN &A, const RMatXN &C, const RMatXN &F, const RMatXN &L, const RMatXN &N, 
    const RMatXN4 &X, Attenuation3D *att):
Elastic3D(att), mA(A), mC(C), mF(F), mL(L), mN(N), 
mN2(two * mN), mX(X) {
    int Nr = mA.rows();
    int n = Nr / 2 + 1;
    mSin1t = (theta.array().sin().matrix()).cast<Real>().colwise().replicate(n);
    mCos1t = (theta.array().cos().matrix()).cast<Real>().colwise().replicate(n);
    mSin2t = ((2. * theta).array().sin().matrix()).cast<Real>().colwise().replicate(n);
    mCos2t = ((2. * theta).array().cos().matrix()).cast<Real>().colwise().replicate(n);
}

void TransverselyIsotropic3X::strainToStress(const vec_ar9_CMatPP &strain, vec_ar9_CMatPP &stress, int dummy) const {
    // constants
    int Nr =  mA.rows();
    int Nu = Nr / 2;
    
    // copy forward
    CMatXN9 &strainCopy = SolverFFTW_N9::getR2C_CMat(Nr);
    Elastic3D::flattenVector(strain, strainCopy, Nu);
    
    // rotate forward
    CMatXN9 &strainTIsoC = SolverFFTW_N9::getC2R_CMat(Nr);
    rotateStrainToTIso(strainCopy, strainTIsoC, Nu);
    
    // FFT forward
    SolverFFTW_N9::computeC2R(Nr);
    RMatXN9 &strainTIsoR9 = SolverFFTW_N9::getC2R_RMat(Nr);
    
    // particle relabelling forward
    RMatXN6 &strainTIsoR6 = SolverFFTW_N6::getC2R_RMat(Nr);
    transformStrainX(strainTIsoR9, strainTIsoR6, Nr);
    
    // strain => stress
    RMatXN6 &stressTIsoR6 = SolverFFTW_N6::getR2C_RMat(Nr);
    
    // to avoid dynamic allocation, use stressTIsoR6.block(0, 3&4 * nPE, Nr, nPE) as temp memory
    stressTIsoR6.block(0, 3 * nPE, Nr, nPE) = strainTIsoR6.block(0, 0 * nPE, Nr, nPE) + strainTIsoR6.block(0, 1 * nPE, Nr, nPE);
    stressTIsoR6.block(0, 4 * nPE, Nr, nPE) = mA.schur(stressTIsoR6.block(0, 3 * nPE, Nr, nPE)) + mF.schur(strainTIsoR6.block(0, 2 * nPE, Nr, nPE));
    
    stressTIsoR6.block(0, 0 * nPE, Nr, nPE) = stressTIsoR6.block(0, 4 * nPE, Nr, nPE) - mN2.schur(strainTIsoR6.block(0, 1 * nPE, Nr, nPE));
    stressTIsoR6.block(0, 1 * nPE, Nr, nPE) = stressTIsoR6.block(0, 4 * nPE, Nr, nPE) - mN2.schur(strainTIsoR6.block(0, 0 * nPE, Nr, nPE));
    stressTIsoR6.block(0, 2 * nPE, Nr, nPE) = mC.schur(strainTIsoR6.block(0, 2 * nPE, Nr, nPE)) + mF.schur(stressTIsoR6.block(0, 3 * nPE, Nr, nPE));
    stressTIsoR6.block(0, 3 * nPE, Nr, nPE) = mL.schur(strainTIsoR6.block(0, 3 * nPE, Nr, nPE));
    stressTIsoR6.block(0, 4 * nPE, Nr, nPE) = mL.schur(strainTIsoR6.block(0, 4 * nPE, Nr, nPE));
    stressTIsoR6.block(0, 5 * nPE, Nr, nPE) = mN.schur(strainTIsoR6.block(0, 5 * nPE, Nr, nPE));
    
    if (mAttenuation) {
        mAttenuation->applyToStress(stressTIsoR6);
        mAttenuation->updateMemoryVariables(strainTIsoR6);
    }
        
    // particle relabelling backward
    RMatXN9 &stressTIsoR9 = SolverFFTW_N9::getR2C_RMat(Nr);
    transformStressX(stressTIsoR6, stressTIsoR9, Nr);
    
    // FFT backward
    SolverFFTW_N9::computeR2C(Nr);
    CMatXN9 &stressTIsoC = SolverFFTW_N9::getR2C_CMat(Nr);
    
    // rotate backward
    CMatXN9 &stressCopy = SolverFFTW_N9::getC2R_CMat(Nr);
    rotateStressToCyln(stressTIsoC, stressCopy, Nu);
    
    // copy backward
    Elastic3D::stackupVector(stressCopy, stress, Nu);
}

void TransverselyIsotropic3X::transformStrainX(const RMatXN9 &strain9, RMatXN6 &strain6, int n) const {
    strain6.block(0, nPE * 0, n, nPE) = mX.block(0, nPE * 0, n, nPE).schur(strain9.block(0, nPE * 0, n, nPE))
                                      + mX.block(0, nPE * 1, n, nPE).schur(strain9.block(0, nPE * 2, n, nPE));
    strain6.block(0, nPE * 1, n, nPE) = mX.block(0, nPE * 0, n, nPE).schur(strain9.block(0, nPE * 4, n, nPE))
                                      + mX.block(0, nPE * 2, n, nPE).schur(strain9.block(0, nPE * 5, n, nPE));
    strain6.block(0, nPE * 2, n, nPE) = mX.block(0, nPE * 3, n, nPE).schur(strain9.block(0, nPE * 8, n, nPE));
    strain6.block(0, nPE * 3, n, nPE) = mX.block(0, nPE * 0, n, nPE).schur(strain9.block(0, nPE * 7, n, nPE))
                                      + mX.block(0, nPE * 2, n, nPE).schur(strain9.block(0, nPE * 8, n, nPE))
                                      + mX.block(0, nPE * 3, n, nPE).schur(strain9.block(0, nPE * 5, n, nPE));
    strain6.block(0, nPE * 4, n, nPE) = mX.block(0, nPE * 0, n, nPE).schur(strain9.block(0, nPE * 6, n, nPE))
                                      + mX.block(0, nPE * 1, n, nPE).schur(strain9.block(0, nPE * 8, n, nPE))
                                      + mX.block(0, nPE * 3, n, nPE).schur(strain9.block(0, nPE * 2, n, nPE));
    strain6.block(0, nPE * 5, n, nPE) = mX.block(0, nPE * 0, n, nPE).schur(strain9.block(0, nPE * 3, n, nPE) + strain9.block(0, nPE * 1, n, nPE))
                                      + mX.block(0, nPE * 1, n, nPE).schur(strain9.block(0, nPE * 5, n, nPE))
                                      + mX.block(0, nPE * 2, n, nPE).schur(strain9.block(0, nPE * 2, n, nPE));
}

void TransverselyIsotropic3X::transformStressX(const RMatXN6 &stress6, RMatXN9 &stress9, int n) const {
    stress9.block(0, nPE * 0, n, nPE) = mX.block(0, nPE * 0, n, nPE).schur(stress6.block(0, nPE * 0, n, nPE));
    stress9.block(0, nPE * 1, n, nPE) = mX.block(0, nPE * 0, n, nPE).schur(stress6.block(0, nPE * 5, n, nPE));
    stress9.block(0, nPE * 2, n, nPE) = mX.block(0, nPE * 1, n, nPE).schur(stress6.block(0, nPE * 0, n, nPE)) 
                                      + mX.block(0, nPE * 3, n, nPE).schur(stress6.block(0, nPE * 4, n, nPE))
                                      + mX.block(0, nPE * 2, n, nPE).schur(stress6.block(0, nPE * 5, n, nPE));
    stress9.block(0, nPE * 3, n, nPE) = stress9.block(0, nPE * 1, n, nPE);                                  
    stress9.block(0, nPE * 4, n, nPE) = mX.block(0, nPE * 0, n, nPE).schur(stress6.block(0, nPE * 1, n, nPE));
    stress9.block(0, nPE * 5, n, nPE) = mX.block(0, nPE * 2, n, nPE).schur(stress6.block(0, nPE * 1, n, nPE)) 
                                      + mX.block(0, nPE * 3, n, nPE).schur(stress6.block(0, nPE * 3, n, nPE))
                                      + mX.block(0, nPE * 1, n, nPE).schur(stress6.block(0, nPE * 5, n, nPE));
    stress9.block(0, nPE * 6, n, nPE) = mX.block(0, nPE * 0, n, nPE).schur(stress6.block(0, nPE * 4, n, nPE));
    stress9.block(0, nPE * 7, n, nPE) = mX.block(0, nPE * 0, n, nPE).schur(stress6.block(0, nPE * 3, n, nPE));
    stress9.block(0, nPE * 8, n, nPE) = mX.block(0, nPE * 3, n, nPE).schur(stress6.block(0, nPE * 2, n, nPE)) 
                                      + mX.block(0, nPE * 2, n, nPE).schur(stress6.block(0, nPE * 3, n, nPE))
                                      + mX.block(0, nPE * 1, n, nPE).schur(stress6.block(0, nPE * 4, n, nPE));
}

void TransverselyIsotropic3X::rotateStrainToTIso(const CMatXN9 &strainCyln, CMatXN9 &strainTIso, int Nu) const {
    int n = Nu + 1;
    // temp
    strainTIso.block(0, nPE * 1, n, nPE) = strainCyln.block(0, nPE * 0, n, nPE) + strainCyln.block(0, nPE * 8, n, nPE); 
    strainTIso.block(0, nPE * 3, n, nPE) = strainCyln.block(0, nPE * 0, n, nPE) - strainCyln.block(0, nPE * 8, n, nPE); 
    strainTIso.block(0, nPE * 5, n, nPE) = strainCyln.block(0, nPE * 2, n, nPE) + strainCyln.block(0, nPE * 6, n, nPE); 
    strainTIso.block(0, nPE * 7, n, nPE) = strainCyln.block(0, nPE * 2, n, nPE) - strainCyln.block(0, nPE * 6, n, nPE); 
    // second order
    strainTIso.block(0, nPE * 0, n, nPE) = half * (strainTIso.block(0, nPE * 1, n, nPE) + mCos2t.schur(strainTIso.block(0, nPE * 3, n, nPE)) - mSin2t.schur(strainTIso.block(0, nPE * 5, n, nPE)));
    strainTIso.block(0, nPE * 2, n, nPE) = half * (strainTIso.block(0, nPE * 7, n, nPE) + mCos2t.schur(strainTIso.block(0, nPE * 5, n, nPE)) + mSin2t.schur(strainTIso.block(0, nPE * 3, n, nPE)));
    strainTIso.block(0, nPE * 6, n, nPE) = strainTIso.block(0, nPE * 2, n, nPE) - strainTIso.block(0, nPE * 7, n, nPE);
    strainTIso.block(0, nPE * 8, n, nPE) = strainTIso.block(0, nPE * 1, n, nPE) - strainTIso.block(0, nPE * 0, n, nPE);
    // first order
    strainTIso.block(0, nPE * 1, n, nPE) = mCos1t.schur(strainCyln.block(0, nPE * 1, n, nPE)) - mSin1t.schur(strainCyln.block(0, nPE * 7, n, nPE));
    strainTIso.block(0, nPE * 7, n, nPE) = mCos1t.schur(strainCyln.block(0, nPE * 7, n, nPE)) + mSin1t.schur(strainCyln.block(0, nPE * 1, n, nPE));    
    strainTIso.block(0, nPE * 3, n, nPE) = mCos1t.schur(strainCyln.block(0, nPE * 3, n, nPE)) - mSin1t.schur(strainCyln.block(0, nPE * 5, n, nPE));
    strainTIso.block(0, nPE * 5, n, nPE) = mCos1t.schur(strainCyln.block(0, nPE * 5, n, nPE)) + mSin1t.schur(strainCyln.block(0, nPE * 3, n, nPE));                             
    // zeroth order
    strainTIso.block(0, nPE * 4, n, nPE) = strainCyln.block(0, nPE * 4, n, nPE);
}

void TransverselyIsotropic3X::rotateStressToCyln(const CMatXN9 &stressTIso, CMatXN9 &stressCyln, int Nu) const
{
    int n = Nu + 1;
    // temp
    stressCyln.block(0, nPE * 1, n, nPE) = stressTIso.block(0, nPE * 0, n, nPE) + stressTIso.block(0, nPE * 8, n, nPE);
    stressCyln.block(0, nPE * 3, n, nPE) = stressTIso.block(0, nPE * 0, n, nPE) - stressTIso.block(0, nPE * 8, n, nPE);
    stressCyln.block(0, nPE * 5, n, nPE) = stressTIso.block(0, nPE * 2, n, nPE) + stressTIso.block(0, nPE * 6, n, nPE);
    stressCyln.block(0, nPE * 7, n, nPE) = stressTIso.block(0, nPE * 2, n, nPE) - stressTIso.block(0, nPE * 6, n, nPE);
    // second order 
    stressCyln.block(0, nPE * 0, n, nPE) = half * (stressCyln.block(0, nPE * 1, n, nPE) + mCos2t.schur(stressCyln.block(0, nPE * 3, n, nPE)) + mSin2t.schur(stressCyln.block(0, nPE * 5, n, nPE)));
    stressCyln.block(0, nPE * 2, n, nPE) = half * (stressCyln.block(0, nPE * 7, n, nPE) + mCos2t.schur(stressCyln.block(0, nPE * 5, n, nPE)) - mSin2t.schur(stressCyln.block(0, nPE * 3, n, nPE)));
    stressCyln.block(0, nPE * 6, n, nPE) = stressCyln.block(0, nPE * 2, n, nPE) - stressCyln.block(0, nPE * 7, n, nPE);
    stressCyln.block(0, nPE * 8, n, nPE) = stressCyln.block(0, nPE * 1, n, nPE) - stressCyln.block(0, nPE * 0, n, nPE);
    // first order
    stressCyln.block(0, nPE * 1, n, nPE) = mCos1t.schur(stressTIso.block(0, nPE * 1, n, nPE)) + mSin1t.schur(stressTIso.block(0, nPE * 7, n, nPE));
    stressCyln.block(0, nPE * 7, n, nPE) = mCos1t.schur(stressTIso.block(0, nPE * 7, n, nPE)) - mSin1t.schur(stressTIso.block(0, nPE * 1, n, nPE));
    stressCyln.block(0, nPE * 3, n, nPE) = mCos1t.schur(stressTIso.block(0, nPE * 3, n, nPE)) + mSin1t.schur(stressTIso.block(0, nPE * 5, n, nPE));
    stressCyln.block(0, nPE * 5, n, nPE) = mCos1t.schur(stressTIso.block(0, nPE * 5, n, nPE)) - mSin1t.schur(stressTIso.block(0, nPE * 3, n, nPE));
    // zeroth order
    stressCyln.block(0, nPE * 4, n, nPE) = stressTIso.block(0, nPE * 4, n, nPE);
}

void TransverselyIsotropic3X::checkCompatibility(int Nr, bool isVoigt) const {
    Elastic3D::checkCompatibility(Nr, isVoigt);
    int myNr = mA.rows();
    if (Nr != myNr) 
        throw std::runtime_error("TransverselyIsotropic3X::checkCompatibility || Incompatible size.");
    if (isVoigt) throw std::runtime_error("TransverselyIsotropic3X::checkCompatibility || "
        "Incompatible gradient operator.");
}

