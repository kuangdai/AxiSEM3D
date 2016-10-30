// Isotropic1D.cpp
// created by Kuangdai on 3-Apr-2016 
// isotropic 1D material

#include "Isotropic1D.h"
#include "Attenuation1D.h"

Isotropic1D::Isotropic1D(const RMatPP &lambda, const RMatPP &mu, Attenuation1D *att): 
Elastic1D(att), mLambda(lambda), mMu(mu), mMu2(two * mu) {
    // nothing
}

void Isotropic1D::strainToStress(const vec_ar9_CMatPP &strain, vec_ar9_CMatPP &stress, int Nu) const {
    static CMatPP sii;
    for (int alpha = 0; alpha <= Nu; alpha++) {
        sii = mLambda.schur(strain[alpha][0] + strain[alpha][1] + strain[alpha][2]);
        stress[alpha][0] = sii + mMu2.schur(strain[alpha][0]);
        stress[alpha][1] = sii + mMu2.schur(strain[alpha][1]);
        stress[alpha][2] = sii + mMu2.schur(strain[alpha][2]);
        stress[alpha][3] = mMu.schur(strain[alpha][3]);
        stress[alpha][4] = mMu.schur(strain[alpha][4]);
        stress[alpha][5] = mMu.schur(strain[alpha][5]);
    }
    if (mAttenuation) {
        mAttenuation->applyToStress(stress);
        mAttenuation->updateMemoryVariables(strain);
    }
}

