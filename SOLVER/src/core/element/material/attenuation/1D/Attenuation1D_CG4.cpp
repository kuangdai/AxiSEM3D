// Attenuation1D_CG44.h
// created by Kuangdai on 29-Apr-2016 
// 1D attenuation on coarse grid


#include "Attenuation1D_CG4.h"

Attenuation1D_CG4::Attenuation1D_CG4(int nsls, const RColX &alpha, 
    const RColX &beta, const RColX &gamma, int Nu, 
    const RRow4 &dkappa, const RRow4 &dmu, bool doKappa):
Attenuation1D(nsls, alpha, beta, gamma),
mDKappa3(three * dkappa), mDMu(dmu), mDMu2(two * dmu), mDoKappa(doKappa) {
    mStressR = vec_ar6_CRow4(Nu + 1, zero_ar6_CRow4);
    mMemVar = std::vector<vec_ar6_CRow4>(mNSLS, mStressR);
}

void Attenuation1D_CG4::applyToStress(vec_ar6_CMatPP &stress) const {
    int Nu = mStressR.size() - 1;
    for (int isls = 0; isls < mNSLS; isls++) {
        for (int alpha = 0; alpha <= Nu; alpha++) {
            for (int i = 0; i < 6; i++) {
                stress[alpha][i](1, 1) -= mMemVar[isls][alpha][i](0);
                stress[alpha][i](1, 3) -= mMemVar[isls][alpha][i](1);
                stress[alpha][i](3, 1) -= mMemVar[isls][alpha][i](2);
                stress[alpha][i](3, 3) -= mMemVar[isls][alpha][i](3);
            }
        }
    }
}

void Attenuation1D_CG4::updateMemoryVariables(const vec_ar6_CMatPP &strain) {
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
    
    static ar6_CRow4 strain4;
    static CRow4 eii_over_3, sii_over_3;
    for (int alpha = 0; alpha <= Nu; alpha++) {
        for (int i = 0; i < 6; i++) {
            strain4[i](0) = strain[alpha][i](1, 1);
            strain4[i](1) = strain[alpha][i](1, 3);
            strain4[i](2) = strain[alpha][i](3, 1);
            strain4[i](3) = strain[alpha][i](3, 3);
        }
        eii_over_3 = (strain4[0] + strain4[1] + strain4[2]) * third;
        if (mDoKappa) {
            sii_over_3 = mDKappa3.schur(eii_over_3);
            mStressR[alpha][0] = sii_over_3 + mDMu2.schur(strain4[0] - eii_over_3);
            mStressR[alpha][1] = sii_over_3 + mDMu2.schur(strain4[1] - eii_over_3);
            mStressR[alpha][2] = sii_over_3 + mDMu2.schur(strain4[2] - eii_over_3);
        } else {
            mStressR[alpha][0] = mDMu2.schur(strain4[0] - eii_over_3);
            mStressR[alpha][1] = mDMu2.schur(strain4[1] - eii_over_3);
            mStressR[alpha][2] = -(mStressR[alpha][0] + mStressR[alpha][1]);
        }
        mStressR[alpha][3] = mDMu.schur(strain4[3]);
        mStressR[alpha][4] = mDMu.schur(strain4[4]);
        mStressR[alpha][5] = mDMu.schur(strain4[5]);
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

void Attenuation1D_CG4::checkCompatibility(int Nr) const
{
    if (nPol != 4) {
        throw std::runtime_error("Attenuation1D_CG4::checkCompatibility "
            " || Coarse grid attenuation is available only when nPol = 4.");
    }
            
    if (Nr / 2 + 1 != mStressR.size()) {
        throw std::runtime_error("Attenuation1D_CG4::checkCompatibility || Incompatible size.");
    }
}

void Attenuation1D_CG4::resetZero() {
    mStressR = vec_ar6_CRow4(mStressR.size(), zero_ar6_CRow4);
    mMemVar = std::vector<vec_ar6_CRow4>(mNSLS, mStressR);
}


