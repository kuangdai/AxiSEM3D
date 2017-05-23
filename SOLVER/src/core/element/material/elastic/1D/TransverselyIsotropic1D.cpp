// TransverselyIsotropic1D.cpp
// created by Kuangdai on 22-Apr-2016 
// transversely isotropic 1D material

#include "TransverselyIsotropic1D.h"
#include "Attenuation1D.h"
#include "SolidElement.h"

void TransverselyIsotropic1D::strainToStress(SolidResponse &response) const {
    static CMatPP e0_p_e1, temp;
    for (int alpha = 0; alpha <= response.mNu; alpha++) {
        const ar6_CMatPP &strainTIso = response.mStrain6[alpha];
        ar6_CMatPP &stressTIso = response.mStress6[alpha];
        e0_p_e1 = strainTIso[0] + strainTIso[1];
        temp = mA.schur(e0_p_e1) + mF.schur(strainTIso[2]);
        stressTIso[0] = temp - mN2.schur(strainTIso[1]);
        stressTIso[1] = temp - mN2.schur(strainTIso[0]);
        stressTIso[2] = mC.schur(strainTIso[2]) + mF.schur(e0_p_e1);
        stressTIso[3] = mL.schur(strainTIso[3]);
        stressTIso[4] = mL.schur(strainTIso[4]);
        stressTIso[5] = mN.schur(strainTIso[5]);
    }
    if (mAttenuation) {
        mAttenuation->applyToStress(response.mStress6);
        mAttenuation->updateMemoryVariables(response.mStrain6);
    }
}
