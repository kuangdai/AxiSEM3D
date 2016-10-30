// Acoustic3X.cpp
// created by Kuangdai on 6-Jun-2016 
// 3D acoustic with topography

#include "Acoustic3X.h"
#include "SolverFFTW_N3.h"

Acoustic3X::Acoustic3X(const RDRowN &theta, const RMatXN KJ, const RMatXN4 X) {
    int n = KJ.rows();
    mSint = (theta.array().sin().matrix()).cast<Real>().colwise().replicate(n);
    mCost = (theta.array().cos().matrix()).cast<Real>().colwise().replicate(n);
    mX00  = KJ.schur(X.block(0, nPE * 0, n, nPE)).schur(X.block(0, nPE * 0, n, nPE));
    mX01  = KJ.schur(X.block(0, nPE * 0, n, nPE)).schur(X.block(0, nPE * 1, n, nPE));
    mX02  = KJ.schur(X.block(0, nPE * 0, n, nPE)).schur(X.block(0, nPE * 2, n, nPE));
    mX123 = KJ.schur(X.block(0, nPE * 1, n, nPE)).schur(X.block(0, nPE * 1, n, nPE)) + 
            KJ.schur(X.block(0, nPE * 2, n, nPE)).schur(X.block(0, nPE * 2, n, nPE)) +
            KJ.schur(X.block(0, nPE * 3, n, nPE)).schur(X.block(0, nPE * 3, n, nPE));
}

void Acoustic3X::strainToStress(const vec_ar3_CMatPP &strain, vec_ar3_CMatPP &stress, int dummy) const {
    // constants
    int Nr = mX00.rows();
    int Nu = Nr / 2;
    
    // copy forward
    CMatXN3 &strainCopy = SolverFFTW_N3::getR2C_CMat(Nr);
    Acoustic3X::flattenScalar(strain, strainCopy, Nu);
    
    // rotate forward
    CMatXN3 &strainTIsoC = SolverFFTW_N3::getC2R_CMat(Nr);
    rotateStrainToTIso(strainCopy, strainTIsoC, Nu);
    
    // FFT forward
    SolverFFTW_N3::computeC2R(Nr);
    RMatXN3 &strainTIsoR = SolverFFTW_N3::getC2R_RMat(Nr);
    RMatXN3 &stressTIsoR = SolverFFTW_N3::getR2C_RMat(Nr);
    
    // strain => stress
    stressTIsoR.block(0, 0 * nPE, Nr, nPE) = mX00.schur(strainTIsoR.block(0, 0 * nPE, Nr, nPE)) +
                                             mX01.schur(strainTIsoR.block(0, 2 * nPE, Nr, nPE));
    stressTIsoR.block(0, 1 * nPE, Nr, nPE) = mX00.schur(strainTIsoR.block(0, 1 * nPE, Nr, nPE)) +
                                             mX02.schur(strainTIsoR.block(0, 2 * nPE, Nr, nPE));
    stressTIsoR.block(0, 2 * nPE, Nr, nPE) = mX01.schur(strainTIsoR.block(0, 0 * nPE, Nr, nPE)) +  
                                             mX02.schur(strainTIsoR.block(0, 1 * nPE, Nr, nPE)) +
                                            mX123.schur(strainTIsoR.block(0, 2 * nPE, Nr, nPE));
    
    // FFT backward
    SolverFFTW_N3::computeR2C(Nr);
    CMatXN3 &stressTIsoC = SolverFFTW_N3::getR2C_CMat(Nr);
    
    // rotate backward
    CMatXN3 &stressCopy = SolverFFTW_N3::getC2R_CMat(Nr);
    rotateStressToCyln(stressTIsoC, stressCopy, Nu);
    
    // copy backward
    Acoustic3X::stackupScalar(stressCopy, stress, Nu);
}

void Acoustic3X::rotateStrainToTIso(const CMatXN3 &strainCyln, CMatXN3 &strainTIso, int Nu) const {
    int n = Nu + 1;
    strainTIso.block(0, nPE * 0, n, nPE) = mCost.schur(strainCyln.block(0, nPE * 0, n, nPE)) 
                                         - mSint.schur(strainCyln.block(0, nPE * 2, n, nPE));
    strainTIso.block(0, nPE * 1, n, nPE) = strainCyln.block(0, nPE * 1, n, nPE);
    strainTIso.block(0, nPE * 2, n, nPE) = mCost.schur(strainCyln.block(0, nPE * 2, n, nPE)) 
                                         + mSint.schur(strainCyln.block(0, nPE * 0, n, nPE));
}

void Acoustic3X::rotateStressToCyln(const CMatXN3 &stressTIso, CMatXN3 &stressCyln, int Nu) const {
    int n = Nu + 1;
    stressCyln.block(0, nPE * 0, n, nPE) = mCost.schur(stressTIso.block(0, nPE * 0, n, nPE)) 
                                         + mSint.schur(stressTIso.block(0, nPE * 2, n, nPE));
    stressCyln.block(0, nPE * 1, n, nPE) = stressTIso.block(0, nPE * 1, n, nPE);
    stressCyln.block(0, nPE * 2, n, nPE) = mCost.schur(stressTIso.block(0, nPE * 2, n, nPE)) 
                                         - mSint.schur(stressTIso.block(0, nPE * 0, n, nPE));
}

void Acoustic3X::checkCompatibility(int Nr) const {
    int myNr = mX00.rows();
    if (Nr != myNr) throw std::runtime_error("Acoustic3X::checkCompatibility || Incompatible size.");
}

void Acoustic3X::flattenScalar(const vec_ar3_CMatPP &mat, CMatXN3 &row, int Nu) {
    for (int alpha = 0; alpha <= Nu; alpha++) 
        for (int i = 0; i < 3; i++) 
            for (int j = 0; j < nPntEdge; j++)
                row.block(alpha, nPE * i + nPntEdge * j, 1, nPntEdge) 
                    = mat[alpha][i].block(j, 0, 1, nPntEdge);
}

void Acoustic3X::stackupScalar(const CMatXN3 &row, vec_ar3_CMatPP &mat, int Nu) {
    for (int alpha = 0; alpha <= Nu; alpha++) 
        for (int i = 0; i < 3; i++) 
            for (int j = 0; j < nPntEdge; j++)
                mat[alpha][i].block(j, 0, 1, nPntEdge) 
                    = row.block(alpha, nPE * i + nPntEdge * j, 1, nPntEdge);
}


