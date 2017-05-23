// Attenuation1D_Full.h
// created by Kuangdai on 29-Apr-2016 
// 1D attenuation on full grid

#include "Attenuation1D_Full.h"

Attenuation1D_Full::Attenuation1D_Full(int nsls, const RColX &alpha, 
    const RColX &beta, const RColX &gamma, int Nu, 
    const RMatPP &dkappa, const RMatPP &dmu, bool doKappa):
Attenuation1D(nsls, alpha, beta, gamma), 
mDKappa3(three * dkappa), mDMu(dmu), mDMu2(two * dmu), mDoKappa(doKappa) {
    mStressR = vec_ar6_CMatPP(Nu + 1, zero_ar6_CMatPP);
    mMemVar = std::vector<vec_ar6_CMatPP>(mNSLS, mStressR);
}

void Attenuation1D_Full::applyToStress(vec_ar6_CMatPP &stress) const {
    int Nu = mStressR.size() - 1;
    for (int isls = 0; isls < mNSLS; isls++) {
        for (int alpha = 0; alpha <= Nu; alpha++) {
            stress[alpha][0] -= mMemVar[isls][alpha][0];
            stress[alpha][1] -= mMemVar[isls][alpha][1];
            stress[alpha][2] -= mMemVar[isls][alpha][2];
            stress[alpha][3] -= mMemVar[isls][alpha][3];
            stress[alpha][4] -= mMemVar[isls][alpha][4];
            stress[alpha][5] -= mMemVar[isls][alpha][5];
        }
    }
}

void Attenuation1D_Full::updateMemoryVariables(const vec_ar6_CMatPP &strain) {
    int Nu = mStressR.size() - 1;
    for (int isls = 0; isls < mNSLS; isls++) {
        Real a = mAlpha[isls];
        Real b = mBeta[isls];
        for (int alpha = 0; alpha <= Nu; alpha++) {
            mMemVar[isls][alpha][0] = a * mMemVar[isls][alpha][0] + b * mStressR[alpha][0];
            mMemVar[isls][alpha][1] = a * mMemVar[isls][alpha][1] + b * mStressR[alpha][1];
            mMemVar[isls][alpha][2] = a * mMemVar[isls][alpha][2] + b * mStressR[alpha][2];
            mMemVar[isls][alpha][3] = a * mMemVar[isls][alpha][3] + b * mStressR[alpha][3];
            mMemVar[isls][alpha][4] = a * mMemVar[isls][alpha][4] + b * mStressR[alpha][4];
            mMemVar[isls][alpha][5] = a * mMemVar[isls][alpha][5] + b * mStressR[alpha][5];
        }
    }
    
    static CMatPP eii_over_3, sii_over_3;
    for (int alpha = 0; alpha <= Nu; alpha++) {
        eii_over_3 = (strain[alpha][0] + strain[alpha][1] + strain[alpha][2]) * third;
        if (mDoKappa) {
            sii_over_3 = mDKappa3.schur(eii_over_3);
            mStressR[alpha][0] = sii_over_3 + mDMu2.schur(strain[alpha][0] - eii_over_3);
            mStressR[alpha][1] = sii_over_3 + mDMu2.schur(strain[alpha][1] - eii_over_3);
            mStressR[alpha][2] = sii_over_3 + mDMu2.schur(strain[alpha][2] - eii_over_3);
        } else {
            mStressR[alpha][0] = mDMu2.schur(strain[alpha][0] - eii_over_3);
            mStressR[alpha][1] = mDMu2.schur(strain[alpha][1] - eii_over_3);
            mStressR[alpha][2] = -(mStressR[alpha][0] + mStressR[alpha][1]);
        }
        mStressR[alpha][3] = mDMu.schur(strain[alpha][3]);
        mStressR[alpha][4] = mDMu.schur(strain[alpha][4]);
        mStressR[alpha][5] = mDMu.schur(strain[alpha][5]);
    }
    
    for (int isls = 0; isls < mNSLS; isls++) {
        Real r = mGamma[isls];
        for (int alpha = 0; alpha <= Nu; alpha++) {
            mMemVar[isls][alpha][0] += r * mStressR[alpha][0];
            mMemVar[isls][alpha][1] += r * mStressR[alpha][1];
            mMemVar[isls][alpha][2] += r * mStressR[alpha][2];
            mMemVar[isls][alpha][3] += r * mStressR[alpha][3];
            mMemVar[isls][alpha][4] += r * mStressR[alpha][4];
            mMemVar[isls][alpha][5] += r * mStressR[alpha][5];
        }
    }
}

void Attenuation1D_Full::checkCompatibility(int Nr) const
{
    if (Nr / 2 + 1 != mStressR.size()) {
        throw std::runtime_error("Attenuation1D_Full::checkCompatibility || Incompatible size.");
    }
}

void Attenuation1D_Full::resetZero() {
    mStressR = vec_ar6_CMatPP(mStressR.size(), zero_ar6_CMatPP);
    mMemVar = std::vector<vec_ar6_CMatPP>(mNSLS, mStressR);
}
