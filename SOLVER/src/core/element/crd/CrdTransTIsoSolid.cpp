// CrdTransTIsoSolid.cpp
// created by Kuangdai on 19-May-2017 
// coordinate transformation between (s, phi, z) and (theta, phi, r)

#include "CrdTransTIsoSolid.h"

CrdTransTIsoSolid::CrdTransTIsoSolid(const RDMatPP &theta) {
    mSin1t = (theta.array().sin().matrix()).cast<Real>();
    mCos1t = (theta.array().cos().matrix()).cast<Real>();
    mSin2t = ((2. * theta).array().sin().matrix()).cast<Real>();
    mCos2t = ((2. * theta).array().cos().matrix()).cast<Real>();
}

void CrdTransTIsoSolid::transformSPZ_RTZ(vec_ar6_CMatPP &u, int Nu) const {
    static CMatPP sum02, dif02, ua_3_;
    for (int alpha = 0; alpha <= Nu; alpha++) {
        ar6_CMatPP &ua = u[alpha];
        sum02 = ua[0] + ua[2];
        dif02 = ua[0] - ua[2];
        ua_3_ = ua[3];
        ua[0] = half * (sum02 + mCos2t.schur(dif02) - mSin2t.schur(ua[4]));
        ua[2] = sum02 - ua[0];
        ua[4] = mCos2t.schur(ua[4]) + mSin2t.schur(dif02);
        ua[3] = mCos1t.schur(ua_3_) + mSin1t.schur(ua[5]);
        ua[5] = mCos1t.schur(ua[5]) - mSin1t.schur(ua_3_);
    }
}

void CrdTransTIsoSolid::transformRTZ_SPZ(vec_ar6_CMatPP &u, int Nu) const {
    static CMatPP sum02, dif02, ua_3_;
    for (int alpha = 0; alpha <= Nu; alpha++) {
        ar6_CMatPP &ua = u[alpha];
        sum02 = ua[0] + ua[2];
        dif02 = (ua[0] - ua[2]) * half;
        ua_3_ = ua[3];
        ua[0] = half * sum02 + mCos2t.schur(dif02) + mSin2t.schur(ua[4]);
        ua[2] = sum02 - ua[0];
        ua[4] = mCos2t.schur(ua[4]) - mSin2t.schur(dif02);
        ua[3] = mCos1t.schur(ua_3_) - mSin1t.schur(ua[5]);
        ua[5] = mCos1t.schur(ua[5]) + mSin1t.schur(ua_3_);
    }
}

void CrdTransTIsoSolid::transformSPZ_RTZ(vec_ar9_CMatPP &u, int Nu) const {
    static CMatPP sum08, dif08, sum26, dif26, ua_1_, ua_3_;
    for (int alpha = 0; alpha <= Nu; alpha++) {
        ar9_CMatPP &ua = u[alpha];
        sum08 = ua[0] + ua[8];
        dif08 = ua[0] - ua[8];
        sum26 = ua[2] + ua[6];
        dif26 = ua[2] - ua[6];
        ua_1_ = ua[1];
        ua_3_ = ua[3];
        ua[0] = half * (sum08 + mCos2t.schur(dif08) - mSin2t.schur(sum26));
        ua[2] = half * (dif26 + mCos2t.schur(sum26) + mSin2t.schur(dif08));
        ua[6] = ua[2] - dif26;
        ua[8] = sum08 - ua[0];
        ua[1] = mCos1t.schur(ua_1_) - mSin1t.schur(ua[7]);
        ua[7] = mCos1t.schur(ua[7]) + mSin1t.schur(ua_1_);    
        ua[3] = mCos1t.schur(ua_3_) - mSin1t.schur(ua[5]);
        ua[5] = mCos1t.schur(ua[5]) + mSin1t.schur(ua_3_);  
    }
}

void CrdTransTIsoSolid::transformRTZ_SPZ(vec_ar9_CMatPP &u, int Nu) const {
    static CMatPP sum08, dif08, sum26, dif26, ua_1_, ua_3_;
    for (int alpha = 0; alpha <= Nu; alpha++) {
        ar9_CMatPP &ua = u[alpha];
        sum08 = ua[0] + ua[8];
        dif08 = ua[0] - ua[8];
        sum26 = ua[2] + ua[6];
        dif26 = ua[2] - ua[6];
        ua_1_ = ua[1];
        ua_3_ = ua[3];
        ua[0] = half * (sum08 + mCos2t.schur(dif08) + mSin2t.schur(sum26));
        ua[2] = half * (dif26 + mCos2t.schur(sum26) - mSin2t.schur(dif08));
        ua[6] = ua[2] - dif26;
        ua[8] = sum08 - ua[0];
        ua[1] = mCos1t.schur(ua_1_) + mSin1t.schur(ua[7]);
        ua[7] = mCos1t.schur(ua[7]) - mSin1t.schur(ua_1_);    
        ua[3] = mCos1t.schur(ua_3_) + mSin1t.schur(ua[5]);
        ua[5] = mCos1t.schur(ua[5]) - mSin1t.schur(ua_3_);  
    }
}

