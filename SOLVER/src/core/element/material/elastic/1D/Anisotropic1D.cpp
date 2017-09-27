// Anisotropic1D.cpp
// created by Kuangdai on 25-Sep-2017 
// full anisotropic 1D material

#include "Anisotropic1D.h"
#include "Attenuation1D.h"
#include "SolidElement.h"

void Anisotropic1D::strainToStress(SolidResponse &response) const {
    for (int alpha = 0; alpha <= response.mNu; alpha++) {
        const ar6_CMatPP &strainTIso = response.mStrain6[alpha];
        ar6_CMatPP &stressTIso = response.mStress6[alpha];
        stressTIso[0] = mC11.schur(strainTIso[0])
                      + mC12.schur(strainTIso[1])
                      + mC13.schur(strainTIso[2])
                      + mC14.schur(strainTIso[3])
                      + mC15.schur(strainTIso[4])
                      + mC16.schur(strainTIso[5]);
        stressTIso[1] = mC12.schur(strainTIso[0])
                      + mC22.schur(strainTIso[1])
                      + mC23.schur(strainTIso[2])
                      + mC24.schur(strainTIso[3])
                      + mC25.schur(strainTIso[4])
                      + mC26.schur(strainTIso[5]);    
        stressTIso[2] = mC13.schur(strainTIso[0])
                      + mC23.schur(strainTIso[1])
                      + mC33.schur(strainTIso[2])
                      + mC34.schur(strainTIso[3])
                      + mC35.schur(strainTIso[4])
                      + mC36.schur(strainTIso[5]);
        stressTIso[3] = mC14.schur(strainTIso[0])
                      + mC24.schur(strainTIso[1])
                      + mC34.schur(strainTIso[2])
                      + mC44.schur(strainTIso[3])
                      + mC45.schur(strainTIso[4])
                      + mC46.schur(strainTIso[5]);
        stressTIso[4] = mC15.schur(strainTIso[0])
                      + mC25.schur(strainTIso[1])
                      + mC35.schur(strainTIso[2])
                      + mC45.schur(strainTIso[3])
                      + mC55.schur(strainTIso[4])
                      + mC56.schur(strainTIso[5]);  
        stressTIso[5] = mC16.schur(strainTIso[0])
                      + mC26.schur(strainTIso[1])
                      + mC36.schur(strainTIso[2])
                      + mC46.schur(strainTIso[3])
                      + mC56.schur(strainTIso[4])
                      + mC66.schur(strainTIso[5]); 
    }
    if (mAttenuation) {
        mAttenuation->applyToStress(response.mStress6);
        mAttenuation->updateMemoryVariables(response.mStrain6);
    }
}
