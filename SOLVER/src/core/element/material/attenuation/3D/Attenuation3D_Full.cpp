// Attenuation3D_Full.h
// created by Kuangdai on 29-Apr-2016 
// 3D attenuation on full grid 

#include "Attenuation3D_Full.h"

Attenuation3D_Full::Attenuation3D_Full(int nsls, 
    const RColX &alpha, const RColX &beta, const RColX &gamma, 
    const RMatXN &dkappa, const RMatXN &dmu, bool doKappa):
Attenuation3D(nsls, alpha, beta, gamma), 
mDKappa3(three * dkappa), mDMu(dmu), mDMu2(two * dmu), mDoKappa(doKappa) {
    mStressR = RMatXN6::Zero(mDMu.rows(), nPE * 6);
    mMemVar = std::vector<RMatXN6>(mNSLS, mStressR);    
}

void Attenuation3D_Full::applyToStress(RMatXN6 &stress) const {
    for (int isls = 0; isls < mNSLS; isls++) {
        stress -= mMemVar[isls];
    } 
}

void Attenuation3D_Full::updateMemoryVariables(const RMatXN6 &strain) {
    for (int isls = 0; isls < mNSLS; isls++) {
        mMemVar[isls] = mAlpha(isls) * mMemVar[isls] + mBeta(isls) * mStressR;
    }
    
    int n = mStressR.rows();
    // to avoid dynamic allocation, use mStressR.block(0, nPE * 3, n, nPE) to store Eii / 3
    mStressR.block(0, nPE * 3, n, nPE) = third * (strain.block(0, nPE * 0, n, nPE) 
        + strain.block(0, nPE * 1, n, nPE) + strain.block(0, nPE * 2, n, nPE));
    if (mDoKappa) {
        // to avoid dynamic allocation, use mStressR.block(0, nPE * 4, n, nPE) to store Sii 
        mStressR.block(0, nPE * 4, n, nPE) = mDKappa3.schur(mStressR.block(0, nPE * 3, n, nPE));
        mStressR.block(0, nPE * 0, n, nPE) = mStressR.block(0, nPE * 4, n, nPE) + mDMu2.schur(strain.block(0, nPE * 0, n, nPE) - mStressR.block(0, nPE * 3, n, nPE));
        mStressR.block(0, nPE * 1, n, nPE) = mStressR.block(0, nPE * 4, n, nPE) + mDMu2.schur(strain.block(0, nPE * 1, n, nPE) - mStressR.block(0, nPE * 3, n, nPE));
        mStressR.block(0, nPE * 2, n, nPE) = mStressR.block(0, nPE * 4, n, nPE) + mDMu2.schur(strain.block(0, nPE * 2, n, nPE) - mStressR.block(0, nPE * 3, n, nPE));
    } else {
        mStressR.block(0, nPE * 0, n, nPE) = mDMu2.schur(strain.block(0, nPE * 0, n, nPE) - mStressR.block(0, nPE * 3, n, nPE));
        mStressR.block(0, nPE * 1, n, nPE) = mDMu2.schur(strain.block(0, nPE * 1, n, nPE) - mStressR.block(0, nPE * 3, n, nPE));
        mStressR.block(0, nPE * 2, n, nPE) = -(mStressR.block(0, nPE * 0, n, nPE) + mStressR.block(0, nPE * 1, n, nPE));
    }
    mStressR.block(0, nPE * 3, n, nPE) = mDMu.schur(strain.block(0, nPE * 3, n, nPE));
    mStressR.block(0, nPE * 4, n, nPE) = mDMu.schur(strain.block(0, nPE * 4, n, nPE));
    mStressR.block(0, nPE * 5, n, nPE) = mDMu.schur(strain.block(0, nPE * 5, n, nPE));
    
    for (int isls = 0; isls < mNSLS; isls++) {
        mMemVar[isls] += mGamma[isls] * mStressR;        
    }
}

void Attenuation3D_Full::checkCompatibility(int Nr) const {
    int myNr = mStressR.rows();
    if (Nr != myNr) {
        throw std::runtime_error("Attenuation3D_Full::checkCompatibility || Incompatible size.");
    }
}

void Attenuation3D_Full::resetZero() {
    mStressR.setZero();
    mMemVar = std::vector<RMatXN6>(mNSLS, mStressR);    
}

