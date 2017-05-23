// Acoustic1D.cpp
// created by Kuangdai on 23-Apr-2016 
// 1D acoustic

#include "Acoustic1D.h"
#include "FluidElement.h"

void Acoustic1D::strainToStress(FluidResponse &response) const {
    for (int alpha = 0; alpha <= response.mNu; alpha++) {
        const ar3_CMatPP &strain_alpha = response.mStrain[alpha];
        ar3_CMatPP &stress_alpha = response.mStress[alpha];
        stress_alpha[0] = mKStruct.schur(strain_alpha[0]);
        stress_alpha[1] = mKStruct.schur(strain_alpha[1]);
        stress_alpha[2] = mKStruct.schur(strain_alpha[2]);
    }
}
