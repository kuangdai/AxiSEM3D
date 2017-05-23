// Attenuation3D_CG4.cpp
// created by Kuangdai on 29-Apr-2016 
// 3D attenuation on coarse grid 

#include "Attenuation3D_CG4.h"

Attenuation3D_CG4::Attenuation3D_CG4(int nsls, 
    const RColX &alpha, const RColX &beta, const RColX &gamma, 
    const RMatX4 &dkappa, const RMatX4 &dmu, bool doKappa):
Attenuation3D(nsls, alpha, beta, gamma), 
mDKappa3(three * dkappa), mDMu(dmu), mDMu2(two * dmu), mDoKappa(doKappa) {
    mStressR = mStrain4 = RMatX46::Zero(mDMu.rows(), nCG * 6);
    mMemVar = std::vector<RMatX46>(mNSLS, mStressR);    
}

void Attenuation3D_CG4::applyToStress(RMatXN6 &stress) const {
    for (int isls = 0; isls < mNSLS; isls++) {
        for (int i = 0; i < 6; i++) {
            stress.col(nPE * i + nPntEdge * 1 + 1) -= mMemVar[isls].col(nCG * i + 0);
            stress.col(nPE * i + nPntEdge * 1 + 3) -= mMemVar[isls].col(nCG * i + 1);
            stress.col(nPE * i + nPntEdge * 3 + 1) -= mMemVar[isls].col(nCG * i + 2);
            stress.col(nPE * i + nPntEdge * 3 + 3) -= mMemVar[isls].col(nCG * i + 3);
        }
    }
}

void Attenuation3D_CG4::updateMemoryVariables(const RMatXN6 &strain) {
    for (int isls = 0; isls < mNSLS; isls++) {
        mMemVar[isls] = mAlpha(isls) * mMemVar[isls] + mBeta(isls) * mStressR;
    }
    
    int n = mStressR.rows();
    for (int i = 0; i < 6; i++) {
        mStrain4.col(nCG * i + 0) = strain.col(nPE * i + nPntEdge * 1 + 1);
        mStrain4.col(nCG * i + 1) = strain.col(nPE * i + nPntEdge * 1 + 3);
        mStrain4.col(nCG * i + 2) = strain.col(nPE * i + nPntEdge * 3 + 1);
        mStrain4.col(nCG * i + 3) = strain.col(nPE * i + nPntEdge * 3 + 3);
    }
    
    // to avoid dynamic allocation, use mStressR.block(0, nCG * 3, n, nCG) to store Eii / 3    
    mStressR.block(0, nCG * 3, n, nCG) = third * (mStrain4.block(0, nCG * 0, n, nCG) 
        + mStrain4.block(0, nCG * 1, n, nCG) + mStrain4.block(0, nCG * 2, n, nCG));
    if (mDoKappa) {
        // to avoid dynamic allocation, use mStressR.block(0, nCG * 4, n, nCG) to store Sii
        mStressR.block(0, nCG * 4, n, nCG) = mDKappa3.schur(mStressR.block(0, nCG * 3, n, nCG));
        mStressR.block(0, nCG * 0, n, nCG) = mStressR.block(0, nCG * 4, n, nCG) + mDMu2.schur(mStrain4.block(0, nCG * 0, n, nCG) - mStressR.block(0, nCG * 3, n, nCG));
        mStressR.block(0, nCG * 1, n, nCG) = mStressR.block(0, nCG * 4, n, nCG) + mDMu2.schur(mStrain4.block(0, nCG * 1, n, nCG) - mStressR.block(0, nCG * 3, n, nCG));
        mStressR.block(0, nCG * 2, n, nCG) = mStressR.block(0, nCG * 4, n, nCG) + mDMu2.schur(mStrain4.block(0, nCG * 2, n, nCG) - mStressR.block(0, nCG * 3, n, nCG));
    } else {
        mStressR.block(0, nCG * 0, n, nCG) = mDMu2.schur(mStrain4.block(0, nCG * 0, n, nCG) - mStressR.block(0, nCG * 3, n, nCG));
        mStressR.block(0, nCG * 1, n, nCG) = mDMu2.schur(mStrain4.block(0, nCG * 1, n, nCG) - mStressR.block(0, nCG * 3, n, nCG));
        mStressR.block(0, nCG * 2, n, nCG) = -(mStressR.block(0, nCG * 0, n, nCG) + mStressR.block(0, nCG * 1, n, nCG));
    }
    mStressR.block(0, nCG * 3, n, nCG) = mDMu.schur(mStrain4.block(0, nCG * 3, n, nCG));
    mStressR.block(0, nCG * 4, n, nCG) = mDMu.schur(mStrain4.block(0, nCG * 4, n, nCG));
    mStressR.block(0, nCG * 5, n, nCG) = mDMu.schur(mStrain4.block(0, nCG * 5, n, nCG));
    
    for (int isls = 0; isls < mNSLS; isls++) {
        mMemVar[isls] += mGamma[isls] * mStressR;
    }
}

void Attenuation3D_CG4::checkCompatibility(int Nr) const
{
    if (nPol != 4) {
        throw std::runtime_error("Attenuation3D_CG4::checkCompatibility "
            " || Coarse grid attenuation is available only when nPol = 4.");
    }
        
    int myNr = mStressR.rows();
    if (Nr != myNr) {
        throw std::runtime_error("Attenuation3D_CG4::checkCompatibility || Incompatible size.");
    }
}

void Attenuation3D_CG4::resetZero() {
    mStressR.setZero();
    mMemVar = std::vector<RMatX46>(mNSLS, mStressR);    
}

