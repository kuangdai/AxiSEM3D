// Isotropic1D.cpp
// created by Kuangdai on 3-Apr-2016 
// isotropic 1D material

#include "Isotropic1D.h"
#include "Attenuation1D.h"
#include "SolidElement.h"

void Isotropic1D::strainToStress(SolidResponse &response) const {
    static CMatPP sii;
    for (int alpha = 0; alpha <= response.mNu; alpha++) {
        const ar6_CMatPP &strain = response.mStrain6[alpha];
        ar6_CMatPP &stress = response.mStress6[alpha];
        sii = mLambda.schur(strain[0] + strain[1] + strain[2]);
        stress[0] = sii + mMu2.schur(strain[0]);
        stress[1] = sii + mMu2.schur(strain[1]);
        stress[2] = sii + mMu2.schur(strain[2]);
        stress[3] = mMu.schur(strain[3]);
        stress[4] = mMu.schur(strain[4]);
        stress[5] = mMu.schur(strain[5]);
    }
    if (mAttenuation) {
        mAttenuation->applyToStress(response.mStress6);
        mAttenuation->updateMemoryVariables(response.mStrain6);
    }
}

