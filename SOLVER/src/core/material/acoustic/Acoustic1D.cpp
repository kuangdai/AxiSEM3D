// Acoustic1D.cpp
// created by Kuangdai on 23-Apr-2016 
// 1D acoustic

#include "Acoustic1D.h"

Acoustic1D::Acoustic1D(const RMatPP &KFluid):
mK(KFluid) {
    // nothing
}

void Acoustic1D::strainToStress(const vec_ar3_CMatPP &strain, vec_ar3_CMatPP &stress, int Nu) const {
    for (int alpha = 0; alpha <= Nu; alpha++) {
        stress[alpha][0] = mK.schur(strain[alpha][0]);
        stress[alpha][1] = mK.schur(strain[alpha][1]);
        stress[alpha][2] = mK.schur(strain[alpha][2]);
    }
}