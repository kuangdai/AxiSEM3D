// PreloopGradient.h
// created by Kuangdai on 20-Sep-2016 
// elemental gradient

#include "PreloopGradient.h"

PreloopGradient::PreloopGradient(const RDMatPP &dsdxii, const RDMatPP &dsdeta, 
                                 const RDMatPP &dzdxii, const RDMatPP &dzdeta, 
                                 const RDMatPP &inv_s, bool axial):
mDsDxii(dsdxii), mDsDeta(dsdeta), 
mDzDxii(dzdxii), mDzDeta(dzdeta), 
mInv_s(inv_s), mAxial(axial) {
    // nothing
}

void PreloopGradient::gradScalar(const vec_CDMatPP &u, vec_ar3_CDMatPP &u_i, int Nu, int nyquist) const {
    static CDMatPP GU, UG;
    for (int alpha = 0; alpha <= Nu - nyquist; alpha++) {
        ComplexD iialpha = (double)alpha * iid;
        GU = (mAxial ? sGT_GLJ : sGT_GLL) * u[alpha];  
        UG = u[alpha] * sG_GLL;
        u_i[alpha][0] = mDzDeta.schur(GU) + mDzDxii.schur(UG);
        u_i[alpha][1] = mInv_s.schur(iialpha * u[alpha]); 
        u_i[alpha][2] = mDsDeta.schur(GU) + mDsDxii.schur(UG);
    }    
    
    // axial masking
    if (mAxial) {
        u_i[0][0].row(0).setZero();
        u_i[0][1].row(0).setZero();
        if (Nu >= 1) {
            u_i[1][1].row(0) = iid * u_i[1][0].row(0);
            u_i[1][2].row(0).setZero();
        }
        for (int alpha = 2; alpha <= Nu; alpha++) {
            u_i[alpha][0].row(0).setZero();
            u_i[alpha][1].row(0).setZero();
            u_i[alpha][2].row(0).setZero();
        }
    }
    
    // mask Nyquist
    if (nyquist) {
        u_i[Nu][0].setZero();
        u_i[Nu][1].setZero();
        u_i[Nu][2].setZero();
    }
}

//-------------------------- static --------------------------//
RDMatPP PreloopGradient::sG_GLL;
RDMatPP PreloopGradient::sG_GLJ;
RDMatPP PreloopGradient::sGT_GLL;
RDMatPP PreloopGradient::sGT_GLJ;
void PreloopGradient::setGMat(const RDMatPP &G_GLL, const RDMatPP &G_GLJ) {
    sG_GLL = G_GLL;
    sG_GLJ = G_GLJ;
    sGT_GLL = sG_GLL.transpose();
    sGT_GLJ = sG_GLJ.transpose();
}

