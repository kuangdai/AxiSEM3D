// Gradient.cpp
// created by Kuangdai on 23-Apr-2016 
// elemental gradient

#include "Gradient.h"

Gradient::Gradient(const RMatPP &dsdxii, const RMatPP &dsdeta, 
                   const RMatPP &dzdxii, const RMatPP &dzdeta, const RMatPP &inv_s):
mDsDxii(dsdxii), mDsDeta(dsdeta), 
mDzDxii(dzdxii), mDzDeta(dzdeta), mInv_s(inv_s) {
    // nothing
}

void Gradient::gradScalar(const vec_CMatPP &u, vec_ar3_CMatPP &u_i, int Nu, int nyquist) const {
    // hardcode for alpha = 0
    static RMatPP GUR, UGR;
    GUR = sGT_GLL * u[0].real();  
    UGR = u[0].real() * sG_GLL;
    u_i[0][0].real() = mDzDeta.schur(GUR) + mDzDxii.schur(UGR);
    u_i[0][2].real() = mDsDeta.schur(GUR) + mDsDxii.schur(UGR);
    
    // alpha > 0
    static CMatPP GU, UG;
    for (int alpha = 1; alpha <= Nu - nyquist; alpha++) {
        Complex iialpha = (Real)alpha * ii;
        GU = sGT_GLL * u[alpha];  
        UG = u[alpha] * sG_GLL;
        u_i[alpha][0] = mDzDeta.schur(GU) + mDzDxii.schur(UG);
        u_i[alpha][1] = mInv_s.schur(iialpha * u[alpha]); 
        u_i[alpha][2] = mDsDeta.schur(GU) + mDsDxii.schur(UG);
    }    
    
    // mask Nyquist
    if (nyquist) {
        u_i[Nu][0].setZero();
        u_i[Nu][1].setZero();
        u_i[Nu][2].setZero();
    }
}

void Gradient::quadScalar(const vec_ar3_CMatPP &f_i, vec_CMatPP &f, int Nu, int nyquist) const {
    // hardcode for mbeta = 0
    static RMatPP XR, YR;
    XR = mDzDeta.schur(f_i[0][0].real()) + mDsDeta.schur(f_i[0][2].real());
    YR = mDzDxii.schur(f_i[0][0].real()) + mDsDxii.schur(f_i[0][2].real());
    f[0].real() = sG_GLL * XR + YR * sGT_GLL; 
    
    // mbeta > 0
    static CMatPP X, Y;
    for (int mbeta = 1; mbeta <= Nu - nyquist; mbeta++) {
        Complex iibeta = - (Real)mbeta * ii; 
        X = mDzDeta.schur(f_i[mbeta][0]) + mDsDeta.schur(f_i[mbeta][2]);
        Y = mDzDxii.schur(f_i[mbeta][0]) + mDsDxii.schur(f_i[mbeta][2]);
        f[mbeta] = sG_GLL * X + Y * sGT_GLL + mInv_s.schur(iibeta * f_i[mbeta][1]);
    }
    
    // mask Nyquist
    if (nyquist) {
        f[Nu].setZero();
    }
}

void Gradient::gradVector(const vec_ar3_CMatPP &ui, vec_ar9_CMatPP &ui_j, int Nu, int nyquist) const {
    // hardcode for alpha = 0
    static RMatPP GU0R, GU1R, GU2R, UG0R, UG1R, UG2R;
    GU0R = sGT_GLL * ui[0][0].real();  
    GU1R = sGT_GLL * ui[0][1].real();  
    GU2R = sGT_GLL * ui[0][2].real();  
    UG0R = ui[0][0].real() * sG_GLL;
    UG1R = ui[0][1].real() * sG_GLL;
    UG2R = ui[0][2].real() * sG_GLL;
    ui_j[0][0].real() = mDzDeta.schur(GU0R) + mDzDxii.schur(UG0R);
    ui_j[0][1].real() = -mInv_s.schur(ui[0][1].real());
    ui_j[0][2].real() = mDsDeta.schur(GU0R) + mDsDxii.schur(UG0R);
    ui_j[0][3].real() = mDzDeta.schur(GU1R) + mDzDxii.schur(UG1R);
    ui_j[0][4].real() = mInv_s.schur(ui[0][0].real()); 
    ui_j[0][5].real() = mDsDeta.schur(GU1R) + mDsDxii.schur(UG1R);
    ui_j[0][6].real() = mDzDeta.schur(GU2R) + mDzDxii.schur(UG2R);
    ui_j[0][8].real() = mDsDeta.schur(GU2R) + mDsDxii.schur(UG2R);
    
    // alpha > 0
    static CMatPP GU0, GU1, GU2, UG0, UG1, UG2;
    for (int alpha = 1; alpha <= Nu - nyquist; alpha++) {        
        Complex iialpha = (Real)alpha * ii;
        GU0 = sGT_GLL * ui[alpha][0];  
        GU1 = sGT_GLL * ui[alpha][1];  
        GU2 = sGT_GLL * ui[alpha][2];  
        UG0 = ui[alpha][0] * sG_GLL;
        UG1 = ui[alpha][1] * sG_GLL;
        UG2 = ui[alpha][2] * sG_GLL;
        ui_j[alpha][0] = mDzDeta.schur(GU0) + mDzDxii.schur(UG0);
        ui_j[alpha][1] = mInv_s.schur(iialpha * ui[alpha][0] - ui[alpha][1]);
        ui_j[alpha][2] = mDsDeta.schur(GU0) + mDsDxii.schur(UG0);
        ui_j[alpha][3] = mDzDeta.schur(GU1) + mDzDxii.schur(UG1);
        ui_j[alpha][4] = mInv_s.schur(ui[alpha][0] + iialpha * ui[alpha][1]); 
        ui_j[alpha][5] = mDsDeta.schur(GU1) + mDsDxii.schur(UG1);
        ui_j[alpha][6] = mDzDeta.schur(GU2) + mDzDxii.schur(UG2);
        ui_j[alpha][7] = mInv_s.schur(iialpha * ui[alpha][2]);
        ui_j[alpha][8] = mDsDeta.schur(GU2) + mDsDxii.schur(UG2);
    }    
    
    // mask Nyquist
    if (nyquist) {
        ui_j[Nu][0].setZero();
        ui_j[Nu][1].setZero();
        ui_j[Nu][2].setZero();
        ui_j[Nu][3].setZero();
        ui_j[Nu][4].setZero();
        ui_j[Nu][5].setZero();
        ui_j[Nu][6].setZero();
        ui_j[Nu][7].setZero();
        ui_j[Nu][8].setZero();
    }
}

void Gradient::quadVector(const vec_ar9_CMatPP &fi_j, vec_ar3_CMatPP &fi, int Nu, int nyquist) const{
    // hardcode for mbeta = 0
    static RMatPP X0R, X1R, X2R, Y0R, Y1R, Y2R;
    X0R = mDzDeta.schur(fi_j[0][0].real()) + mDsDeta.schur(fi_j[0][2].real());
    X1R = mDzDeta.schur(fi_j[0][3].real()) + mDsDeta.schur(fi_j[0][5].real());
    X2R = mDzDeta.schur(fi_j[0][6].real()) + mDsDeta.schur(fi_j[0][8].real());
    Y0R = mDzDxii.schur(fi_j[0][0].real()) + mDsDxii.schur(fi_j[0][2].real());
    Y1R = mDzDxii.schur(fi_j[0][3].real()) + mDsDxii.schur(fi_j[0][5].real());
    Y2R = mDzDxii.schur(fi_j[0][6].real()) + mDsDxii.schur(fi_j[0][8].real());
    fi[0][0].real() = sG_GLL * X0R + Y0R * sGT_GLL + mInv_s.schur(fi_j[0][4].real());
    fi[0][1].real() = sG_GLL * X1R + Y1R * sGT_GLL - mInv_s.schur(fi_j[0][1].real());
    fi[0][2].real() = sG_GLL * X2R + Y2R * sGT_GLL;
    
    // mbeta > 0
    static CMatPP X0, X1, X2, Y0, Y1, Y2;
    for (int mbeta = 1; mbeta <= Nu - nyquist; mbeta++) {
        Complex iibeta = - (Real)mbeta * ii; 
        X0 = mDzDeta.schur(fi_j[mbeta][0]) + mDsDeta.schur(fi_j[mbeta][2]);
        X1 = mDzDeta.schur(fi_j[mbeta][3]) + mDsDeta.schur(fi_j[mbeta][5]);
        X2 = mDzDeta.schur(fi_j[mbeta][6]) + mDsDeta.schur(fi_j[mbeta][8]);
        Y0 = mDzDxii.schur(fi_j[mbeta][0]) + mDsDxii.schur(fi_j[mbeta][2]);
        Y1 = mDzDxii.schur(fi_j[mbeta][3]) + mDsDxii.schur(fi_j[mbeta][5]);
        Y2 = mDzDxii.schur(fi_j[mbeta][6]) + mDsDxii.schur(fi_j[mbeta][8]);
        fi[mbeta][0] = sG_GLL * X0 + Y0 * sGT_GLL + mInv_s.schur(fi_j[mbeta][4] + iibeta * fi_j[mbeta][1]);
        fi[mbeta][1] = sG_GLL * X1 + Y1 * sGT_GLL + mInv_s.schur(iibeta * fi_j[mbeta][4] - fi_j[mbeta][1]);
        fi[mbeta][2] = sG_GLL * X2 + Y2 * sGT_GLL + mInv_s.schur(iibeta * fi_j[mbeta][7]);
    }
    
    // mask Nyquist
    if (nyquist) {
        fi[Nu][0].setZero();
        fi[Nu][1].setZero();
        fi[Nu][2].setZero();
    }
}

//-------------------------- static --------------------------//
RMatPP Gradient::sG_GLL;
RMatPP Gradient::sG_GLJ;
RMatPP Gradient::sGT_GLL;
RMatPP Gradient::sGT_GLJ;
void Gradient::setGMat(const RMatPP &G_GLL, const RMatPP &G_GLJ) {
    sG_GLL = G_GLL;
    sG_GLJ = G_GLJ;
    sGT_GLL = sG_GLL.transpose();
    sGT_GLJ = sG_GLJ.transpose();
}

