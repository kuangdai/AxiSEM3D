// TransverselyIsotropic1D.cpp
// created by Kuangdai on 22-Apr-2016 
// transversely isotropic 1D material

#include "TransverselyIsotropic1D.h"
#include "Attenuation1D.h"

TransverselyIsotropic1D::TransverselyIsotropic1D(const RDMatPP &theta, const RMatPP &A, const RMatPP &C, 
    const RMatPP &F, const RMatPP &L, const RMatPP &N, Attenuation1D *att): 
Elastic1D(att) , mA(A), mC(C), mF(F), mL(L), mN(N), mN2(two * N) {
    mSin1t.topLeftCorner(nPntEdge, nPntEdge) = (theta.array().sin().matrix()).cast<Real>();
    mCos1t.topLeftCorner(nPntEdge, nPntEdge) = (theta.array().cos().matrix()).cast<Real>();
    mSin2t.topLeftCorner(nPntEdge, nPntEdge) = ((2. * theta).array().sin().matrix()).cast<Real>();
    mCos2t.topLeftCorner(nPntEdge, nPntEdge) = ((2. * theta).array().cos().matrix()).cast<Real>();
}

void TransverselyIsotropic1D::strainToStress(const vec_ar9_CMatPP &strain, vec_ar9_CMatPP &stress, int Nu) const {
    static ar6_CMatPP strainTIso, stressTIso;
    static CMatPP e0_p_e1, temp;
    for (int alpha = 0; alpha <= Nu; alpha++) {
        rotateStrainToTIso(strain[alpha], strainTIso);
        e0_p_e1 = strainTIso[0] + strainTIso[1];
        temp = mA.schur(e0_p_e1) + mF.schur(strainTIso[2]);
        stressTIso[0] = temp - mN2.schur(strainTIso[1]);
        stressTIso[1] = temp - mN2.schur(strainTIso[0]);
        stressTIso[2] = mC.schur(strainTIso[2]) + mF.schur(e0_p_e1);
        stressTIso[3] = mL.schur(strainTIso[3]);
        stressTIso[4] = mL.schur(strainTIso[4]);
        stressTIso[5] = mN.schur(strainTIso[5]);
        rotateStressToCyln(stressTIso, stress[alpha]);
    }
    if (mAttenuation) {
        mAttenuation->applyToStress(stress);
        mAttenuation->updateMemoryVariables(strain);
    }
}

void TransverselyIsotropic1D::rotateStrainToTIso(const ar9_CMatPP &strainCyln, ar6_CMatPP &strainTIso) const {
    static CMatPP sum, dif;
    sum = strainCyln[0] + strainCyln[2];
    dif = strainCyln[0] - strainCyln[2];
    strainTIso[0] = half * (sum + mCos2t.schur(dif) - mSin2t.schur(strainCyln[4]));
    strainTIso[1] = strainCyln[1];
    strainTIso[2] = sum - strainTIso[0];
    strainTIso[3] = mCos1t.schur(strainCyln[3]) + mSin1t.schur(strainCyln[5]);
    strainTIso[4] = mCos2t.schur(strainCyln[4]) + mSin2t.schur(dif);
    strainTIso[5] = mCos1t.schur(strainCyln[5]) - mSin1t.schur(strainCyln[3]);
}

void TransverselyIsotropic1D::rotateStressToCyln(const ar6_CMatPP &stressTIso, ar9_CMatPP &stressCyln) const {
    static CMatPP sum, dif;
    sum = stressTIso[0] + stressTIso[2];
    dif = (stressTIso[0] - stressTIso[2]) * half;
    stressCyln[0] = half * sum + mCos2t.schur(dif) + mSin2t.schur(stressTIso[4]);
    stressCyln[1] = stressTIso[1];
    stressCyln[2] = sum - stressCyln[0];
    stressCyln[3] = mCos1t.schur(stressTIso[3]) - mSin1t.schur(stressTIso[5]);
    stressCyln[4] = mCos2t.schur(stressTIso[4]) - mSin2t.schur(dif);
    stressCyln[5] = mCos1t.schur(stressTIso[5]) + mSin1t.schur(stressTIso[3]);
}

