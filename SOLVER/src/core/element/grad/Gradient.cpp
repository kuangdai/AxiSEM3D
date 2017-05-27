// Gradient.cpp
// created by Kuangdai on 19-May-2017 
// elemental gradient

#include "Gradient.h"
#include "FluidElement.h"
#include "SolidElement.h"

Gradient::Gradient(const RDMatPP &dsdxii, const RDMatPP &dsdeta, 
                   const RDMatPP &dzdxii, const RDMatPP &dzdeta, 
                   const RDMatPP &inv_s, bool axial):
mDsDxii(dsdxii.cast<Real>()), mDsDeta(dsdeta.cast<Real>()), 
mDzDxii(dzdxii.cast<Real>()), mDzDeta(dzdeta.cast<Real>()), 
mInv_s(inv_s.cast<Real>()), mAxial(axial) {
    if (mAxial) {
        sG_xii = &sG_GLJ;
        sGT_xii = &sGT_GLJ;
    } else {
        sG_xii = &sG_GLL;
        sGT_xii = &sGT_GLL;
    }
    sG_eta = &sG_GLL;
    sGT_eta = &sGT_GLL;
}

void Gradient::computeGrad(const vec_CMatPP &u, vec_ar3_CMatPP &u_i, int Nu, int nyquist) const {
    // hardcode for alpha = 0
    static RMatPP GUR, UGR;
    GUR = (*sGT_xii) * u[0].real();  
    UGR = u[0].real() * (*sG_eta);
    u_i[0][0].real() = mDzDeta.schur(GUR) + mDzDxii.schur(UGR);
    u_i[0][1].real().setZero();
    u_i[0][2].real() = mDsDeta.schur(GUR) + mDsDxii.schur(UGR);
    
    // alpha > 0
    static CMatPP v, GU, UG;
    for (int alpha = 1; alpha <= Nu - nyquist; alpha++) {        
        Complex iialpha = (Real)alpha * ii;
        v = iialpha * u[alpha];
        GU = (*sGT_xii) * u[alpha];  
        UG = u[alpha] * (*sG_eta);
        u_i[alpha][0] = mDzDeta.schur(GU) + mDzDxii.schur(UG);
        u_i[alpha][1] = mInv_s.schur(v); 
        u_i[alpha][2] = mDsDeta.schur(GU) + mDsDxii.schur(UG);
        if (mAxial) {
            u_i[alpha][1].row(0) += mDzDeta.row(0).schur((*sGT_xii).row(0) * v);
        }
    }    
    
    // mask Nyquist
    if (nyquist) {
        u_i[Nu][0].setZero();
        u_i[Nu][1].setZero();
        u_i[Nu][2].setZero();
    }
}

void Gradient::computeQuad(vec_CMatPP &f, const vec_ar3_CMatPP &f_i, int Nu, int nyquist) const {
    // hardcode for mbeta = 0
    static RMatPP XR, YR;
    XR = mDzDeta.schur(f_i[0][0].real()) + mDsDeta.schur(f_i[0][2].real());
    YR = mDzDxii.schur(f_i[0][0].real()) + mDsDxii.schur(f_i[0][2].real());
    f[0].real() = (*sG_xii) * XR + YR * (*sGT_eta); 
    
    // mbeta > 0
    static CMatPP g, X, Y;
    for (int mbeta = 1; mbeta <= Nu - nyquist; mbeta++) {
        Complex iibeta = - (Real)mbeta * ii; 
        g = iibeta * f_i[mbeta][1];
        X = mDzDeta.schur(f_i[mbeta][0]) + mDsDeta.schur(f_i[mbeta][2]);
        Y = mDzDxii.schur(f_i[mbeta][0]) + mDsDxii.schur(f_i[mbeta][2]);
        f[mbeta] = (*sG_xii) * X + Y * (*sGT_eta) + mInv_s.schur(g);
        if (mAxial) {
            f[mbeta] += (*sG_xii).col(0) * mDzDeta.row(0).schur(g.row(0));
        }
    }
    
    // mask Nyquist
    if (nyquist) {
        f[Nu].setZero();
    }
}    

void Gradient::computeGrad9(const vec_ar3_CMatPP &ui, vec_ar9_CMatPP &ui_j, int Nu, int nyquist) const {
    // hardcode for alpha = 0
    static RMatPP GU0R, GU1R, GU2R, UG0R, UG1R, UG2R;
    GU0R = (*sGT_xii) * ui[0][0].real();  
    GU1R = (*sGT_xii) * ui[0][1].real();  
    GU2R = (*sGT_xii) * ui[0][2].real();  
    UG0R = ui[0][0].real() * (*sG_eta);
    UG1R = ui[0][1].real() * (*sG_eta);
    UG2R = ui[0][2].real() * (*sG_eta);
    ui_j[0][0].real() = mDzDeta.schur(GU0R) + mDzDxii.schur(UG0R);
    ui_j[0][1].real() = -mInv_s.schur(ui[0][1].real());
    ui_j[0][2].real() = mDsDeta.schur(GU0R) + mDsDxii.schur(UG0R);
    ui_j[0][3].real() = mDzDeta.schur(GU1R) + mDzDxii.schur(UG1R);
    ui_j[0][4].real() = mInv_s.schur(ui[0][0].real()); 
    ui_j[0][5].real() = mDsDeta.schur(GU1R) + mDsDxii.schur(UG1R);
    ui_j[0][6].real() = mDzDeta.schur(GU2R) + mDzDxii.schur(UG2R);
    ui_j[0][7].real().setZero();
    ui_j[0][8].real() = mDsDeta.schur(GU2R) + mDsDxii.schur(UG2R);
    if (mAxial) {
        ui_j[0][4].row(0).real() += mDzDeta.row(0).schur((*sGT_xii).row(0) * ui[0][0].real());
        ui_j[0][1].row(0).real() -= mDzDeta.row(0).schur((*sGT_xii).row(0) * ui[0][1].real());
    }
    
    // alpha > 0
    static CMatPP v0, v1, v2, GU0, GU1, GU2, UG0, UG1, UG2;
    for (int alpha = 1; alpha <= Nu - nyquist; alpha++) {        
        Complex iialpha = (Real)alpha * ii;
        v0 = ui[alpha][0] + iialpha * ui[alpha][1];
        v1 = iialpha * ui[alpha][0] - ui[alpha][1];
        v2 = iialpha * ui[alpha][2];
        GU0 = (*sGT_xii) * ui[alpha][0];  
        GU1 = (*sGT_xii) * ui[alpha][1];  
        GU2 = (*sGT_xii) * ui[alpha][2];  
        UG0 = ui[alpha][0] * (*sG_eta);
        UG1 = ui[alpha][1] * (*sG_eta);
        UG2 = ui[alpha][2] * (*sG_eta);
        ui_j[alpha][0] = mDzDeta.schur(GU0) + mDzDxii.schur(UG0);
        ui_j[alpha][1] = mInv_s.schur(v1);
        ui_j[alpha][2] = mDsDeta.schur(GU0) + mDsDxii.schur(UG0);
        ui_j[alpha][3] = mDzDeta.schur(GU1) + mDzDxii.schur(UG1);
        ui_j[alpha][4] = mInv_s.schur(v0); 
        ui_j[alpha][5] = mDsDeta.schur(GU1) + mDsDxii.schur(UG1);
        ui_j[alpha][6] = mDzDeta.schur(GU2) + mDzDxii.schur(UG2);
        ui_j[alpha][7] = mInv_s.schur(v2);
        ui_j[alpha][8] = mDsDeta.schur(GU2) + mDsDxii.schur(UG2);
        if (mAxial) {
            ui_j[alpha][4].row(0) += mDzDeta.row(0).schur((*sGT_xii).row(0) * v0);
            ui_j[alpha][1].row(0) += mDzDeta.row(0).schur((*sGT_xii).row(0) * v1);
            ui_j[alpha][7].row(0) += mDzDeta.row(0).schur((*sGT_xii).row(0) * v2);
            if (alpha == 1) {
                ui_j[alpha][4].row(0) += mDzDxii.row(0).schur(v0.row(0) * (*sG_eta));
                ui_j[alpha][1].row(0) += mDzDxii.row(0).schur(v1.row(0) * (*sG_eta));
            }
        }
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

void Gradient::computeQuad9(vec_ar3_CMatPP &fi, const vec_ar9_CMatPP &fi_j, int Nu, int nyquist) const {
    // hardcode for mbeta = 0
    static RMatPP X0R, X1R, X2R, Y0R, Y1R, Y2R; 
    X0R = mDzDeta.schur(fi_j[0][0].real()) + mDsDeta.schur(fi_j[0][2].real());
    X1R = mDzDeta.schur(fi_j[0][3].real()) + mDsDeta.schur(fi_j[0][5].real());
    X2R = mDzDeta.schur(fi_j[0][6].real()) + mDsDeta.schur(fi_j[0][8].real());
    Y0R = mDzDxii.schur(fi_j[0][0].real()) + mDsDxii.schur(fi_j[0][2].real());
    Y1R = mDzDxii.schur(fi_j[0][3].real()) + mDsDxii.schur(fi_j[0][5].real());
    Y2R = mDzDxii.schur(fi_j[0][6].real()) + mDsDxii.schur(fi_j[0][8].real());
    fi[0][0].real() = (*sG_xii) * X0R + Y0R * (*sGT_eta) + mInv_s.schur(fi_j[0][4].real());
    fi[0][1].real() = (*sG_xii) * X1R + Y1R * (*sGT_eta) - mInv_s.schur(fi_j[0][1].real());
    fi[0][2].real() = (*sG_xii) * X2R + Y2R * (*sGT_eta);
    if (mAxial) {
        fi[0][0].real() += (*sG_xii).col(0) * mDzDeta.row(0).schur(fi_j[0][4].real().row(0));
        fi[0][1].real() -= (*sG_xii).col(0) * mDzDeta.row(0).schur(fi_j[0][1].real().row(0));
    }
    
    // mbeta > 0
    static CMatPP g0, g1, g2, X0, X1, X2, Y0, Y1, Y2;
    for (int mbeta = 1; mbeta <= Nu - nyquist; mbeta++) {
        Complex iibeta = - (Real)mbeta * ii; 
        g0 = fi_j[mbeta][4] + iibeta * fi_j[mbeta][1];
        g1 = iibeta * fi_j[mbeta][4] - fi_j[mbeta][1];
        g2 = iibeta * fi_j[mbeta][7];    
        X0 = mDzDeta.schur(fi_j[mbeta][0]) + mDsDeta.schur(fi_j[mbeta][2]);
        X1 = mDzDeta.schur(fi_j[mbeta][3]) + mDsDeta.schur(fi_j[mbeta][5]);
        X2 = mDzDeta.schur(fi_j[mbeta][6]) + mDsDeta.schur(fi_j[mbeta][8]);
        Y0 = mDzDxii.schur(fi_j[mbeta][0]) + mDsDxii.schur(fi_j[mbeta][2]);
        Y1 = mDzDxii.schur(fi_j[mbeta][3]) + mDsDxii.schur(fi_j[mbeta][5]);
        Y2 = mDzDxii.schur(fi_j[mbeta][6]) + mDsDxii.schur(fi_j[mbeta][8]);
        fi[mbeta][0] = (*sG_xii) * X0 + Y0 * (*sGT_eta) + mInv_s.schur(g0);
        fi[mbeta][1] = (*sG_xii) * X1 + Y1 * (*sGT_eta) + mInv_s.schur(g1);
        fi[mbeta][2] = (*sG_xii) * X2 + Y2 * (*sGT_eta) + mInv_s.schur(g2);
        if (mAxial) {
            fi[mbeta][0] += (*sG_xii).col(0) * mDzDeta.row(0).schur(g0.row(0));
            fi[mbeta][1] += (*sG_xii).col(0) * mDzDeta.row(0).schur(g1.row(0));
            fi[mbeta][2] += (*sG_xii).col(0) * mDzDeta.row(0).schur(g2.row(0));
            if (mbeta == 1) {
                fi[mbeta][0].row(0) += mDzDxii.row(0).schur(g0.row(0)) * (*sGT_eta);
                fi[mbeta][1].row(0) += mDzDxii.row(0).schur(g1.row(0)) * (*sGT_eta);
            }
        }
    }
    
    // mask Nyquist
    if (nyquist) {
        fi[Nu][0].setZero();
        fi[Nu][1].setZero();
        fi[Nu][2].setZero();
    }
}

void Gradient::computeGrad6(const vec_ar3_CMatPP &ui, vec_ar6_CMatPP &eij, int Nu, int nyquist) const {
    // hardcode for alpha = 0
    static RMatPP GU0R, GU1R, GU2R, UG0R, UG1R, UG2R;
    GU0R = (*sGT_xii) * ui[0][0].real();  
    GU1R = (*sGT_xii) * ui[0][1].real();  
    GU2R = (*sGT_xii) * ui[0][2].real();  
    UG0R = ui[0][0].real() * (*sG_eta);
    UG1R = ui[0][1].real() * (*sG_eta);
    UG2R = ui[0][2].real() * (*sG_eta);
    eij[0][0].real() = mDzDeta.schur(GU0R) + mDzDxii.schur(UG0R);
    eij[0][1].real() = mInv_s.schur(ui[0][0].real()); 
    eij[0][2].real() = mDsDeta.schur(GU2R) + mDsDxii.schur(UG2R);
    eij[0][3].real() = mDsDeta.schur(GU1R) + mDsDxii.schur(UG1R);
    eij[0][4].real() = mDsDeta.schur(GU0R) + mDsDxii.schur(UG0R) + mDzDeta.schur(GU2R) + mDzDxii.schur(UG2R);
    eij[0][5].real() = mDzDeta.schur(GU1R) + mDzDxii.schur(UG1R) - mInv_s.schur(ui[0][1].real());
    if (mAxial) {
        eij[0][1].row(0).real() += mDzDeta.row(0).schur((*sGT_xii).row(0) * ui[0][0].real());
        eij[0][5].row(0).real() -= mDzDeta.row(0).schur((*sGT_xii).row(0) * ui[0][1].real());
    }
    
    // alpha > 0
    static CMatPP v0, v1, v2, GU0, GU1, GU2, UG0, UG1, UG2;
    for (int alpha = 1; alpha <= Nu - nyquist; alpha++) {        
        Complex iialpha = (Real)alpha * ii;
        v0 = ui[alpha][0] + iialpha * ui[alpha][1];
        v1 = iialpha * ui[alpha][0] - ui[alpha][1];
        v2 = iialpha * ui[alpha][2];
        GU0 = (*sGT_xii) * ui[alpha][0];  
        GU1 = (*sGT_xii) * ui[alpha][1];  
        GU2 = (*sGT_xii) * ui[alpha][2];  
        UG0 = ui[alpha][0] * (*sG_eta);
        UG1 = ui[alpha][1] * (*sG_eta);
        UG2 = ui[alpha][2] * (*sG_eta);
        eij[alpha][0] = mDzDeta.schur(GU0) + mDzDxii.schur(UG0);
        eij[alpha][1] = mInv_s.schur(v0); 
        eij[alpha][2] = mDsDeta.schur(GU2) + mDsDxii.schur(UG2);
        eij[alpha][3] = mDsDeta.schur(GU1) + mDsDxii.schur(UG1) + mInv_s.schur(v2);
        eij[alpha][4] = mDsDeta.schur(GU0) + mDsDxii.schur(UG0) + mDzDeta.schur(GU2) + mDzDxii.schur(UG2);
        eij[alpha][5] = mDzDeta.schur(GU1) + mDzDxii.schur(UG1) + mInv_s.schur(v1);
        if (mAxial) {
            eij[alpha][1].row(0) += mDzDeta.row(0).schur((*sGT_xii).row(0) * v0);
            eij[alpha][5].row(0) += mDzDeta.row(0).schur((*sGT_xii).row(0) * v1);
            eij[alpha][3].row(0) += mDzDeta.row(0).schur((*sGT_xii).row(0) * v2);
            if (alpha == 1) {
                eij[alpha][1].row(0) += mDzDxii.row(0).schur(v0.row(0) * (*sG_eta));
                eij[alpha][5].row(0) += mDzDxii.row(0).schur(v1.row(0) * (*sG_eta));
            }
        }
    }    
    
    // mask Nyquist
    if (nyquist) {
        eij[Nu][0].setZero();
        eij[Nu][1].setZero();
        eij[Nu][2].setZero();
        eij[Nu][3].setZero();
        eij[Nu][4].setZero();
        eij[Nu][5].setZero();
    }   
}

void Gradient::computeQuad6(vec_ar3_CMatPP &fi, const vec_ar6_CMatPP &sij, int Nu, int nyquist) const {
    // hardcode for mbeta = 0
    static RMatPP X0R, X1R, X2R, Y0R, Y1R, Y2R; 
    X0R = mDzDeta.schur(sij[0][0].real()) + mDsDeta.schur(sij[0][4].real());
    X1R = mDzDeta.schur(sij[0][5].real()) + mDsDeta.schur(sij[0][3].real());
    X2R = mDzDeta.schur(sij[0][4].real()) + mDsDeta.schur(sij[0][2].real());
    Y0R = mDzDxii.schur(sij[0][0].real()) + mDsDxii.schur(sij[0][4].real());
    Y1R = mDzDxii.schur(sij[0][5].real()) + mDsDxii.schur(sij[0][3].real());
    Y2R = mDzDxii.schur(sij[0][4].real()) + mDsDxii.schur(sij[0][2].real());
    fi[0][0].real() = (*sG_xii) * X0R + Y0R * (*sGT_eta) + mInv_s.schur(sij[0][1].real());
    fi[0][1].real() = (*sG_xii) * X1R + Y1R * (*sGT_eta) - mInv_s.schur(sij[0][5].real());
    fi[0][2].real() = (*sG_xii) * X2R + Y2R * (*sGT_eta); 
    if (mAxial) {
        fi[0][0].real() += (*sG_xii).col(0) * mDzDeta.row(0).schur(sij[0][1].real().row(0));
        fi[0][1].real() -= (*sG_xii).col(0) * mDzDeta.row(0).schur(sij[0][5].real().row(0));
    }
    
    // mbeta > 0
    static CMatPP g0, g1, g2, X0, X1, X2, Y0, Y1, Y2;
    for (int mbeta = 1; mbeta <= Nu - nyquist; mbeta++) {
        Complex iibeta = - (Real)mbeta * ii; 
        g0 = sij[mbeta][1] + iibeta * sij[mbeta][5];
        g1 = iibeta * sij[mbeta][1] - sij[mbeta][5];
        g2 = iibeta * sij[mbeta][3];    
        X0 = mDzDeta.schur(sij[mbeta][0]) + mDsDeta.schur(sij[mbeta][4]);
        X1 = mDzDeta.schur(sij[mbeta][5]) + mDsDeta.schur(sij[mbeta][3]);
        X2 = mDzDeta.schur(sij[mbeta][4]) + mDsDeta.schur(sij[mbeta][2]);
        Y0 = mDzDxii.schur(sij[mbeta][0]) + mDsDxii.schur(sij[mbeta][4]);
        Y1 = mDzDxii.schur(sij[mbeta][5]) + mDsDxii.schur(sij[mbeta][3]);
        Y2 = mDzDxii.schur(sij[mbeta][4]) + mDsDxii.schur(sij[mbeta][2]);
        fi[mbeta][0] = (*sG_xii) * X0 + Y0 * (*sGT_eta) + mInv_s.schur(g0);
        fi[mbeta][1] = (*sG_xii) * X1 + Y1 * (*sGT_eta) + mInv_s.schur(g1);
        fi[mbeta][2] = (*sG_xii) * X2 + Y2 * (*sGT_eta) + mInv_s.schur(g2);
        if (mAxial) {
            fi[mbeta][0] += (*sG_xii).col(0) * mDzDeta.row(0).schur(g0.row(0));
            fi[mbeta][1] += (*sG_xii).col(0) * mDzDeta.row(0).schur(g1.row(0));
            fi[mbeta][2] += (*sG_xii).col(0) * mDzDeta.row(0).schur(g2.row(0));
            if (mbeta == 1) {
                fi[mbeta][0].row(0) += mDzDxii.row(0).schur(g0.row(0)) * (*sGT_eta);
                fi[mbeta][1].row(0) += mDzDxii.row(0).schur(g1.row(0)) * (*sGT_eta);
            }
        }
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
void Gradient::setGMat(const RDMatPP &G_GLL, const RDMatPP &G_GLJ) {
    sG_GLL = G_GLL.cast<Real>();
    sG_GLJ = G_GLJ.cast<Real>();
    sGT_GLL = sG_GLL.transpose();
    sGT_GLJ = sG_GLJ.transpose();
}
