// CrdTransTIsoFluid.cpp
// created by Kuangdai on 19-May-2017 
// coordinate transformation between (s, phi, z) and (theta, phi, r)

#include "CrdTransTIsoFluid.h"

CrdTransTIsoFluid::CrdTransTIsoFluid(const RDMatPP &theta) {
    mSin1t = (theta.array().sin().matrix()).cast<Real>();
    mCos1t = (theta.array().cos().matrix()).cast<Real>();
}

void CrdTransTIsoFluid::transformSPZ_RTZ(vec_ar3_CMatPP &u, int Nu) const {
    static CMatPP ua_0_;
    for (int alpha = 0; alpha <= Nu; alpha++) {
        ar3_CMatPP &ua = u[alpha];
        ua_0_ = ua[0];
        ua[0] = mCos1t.schur(ua_0_) - mSin1t.schur(ua[2]);
        ua[2] = mCos1t.schur(ua[2]) + mSin1t.schur(ua_0_);
    }
}

void CrdTransTIsoFluid::transformRTZ_SPZ(vec_ar3_CMatPP &u, int Nu) const {
    static CMatPP ua_0_;
    for (int alpha = 0; alpha <= Nu; alpha++) {
        ar3_CMatPP &ua = u[alpha];
        ua_0_ = ua[0];
        ua[0] = mCos1t.schur(ua_0_) + mSin1t.schur(ua[2]);
        ua[2] = mCos1t.schur(ua[2]) - mSin1t.schur(ua_0_);
    }
}

