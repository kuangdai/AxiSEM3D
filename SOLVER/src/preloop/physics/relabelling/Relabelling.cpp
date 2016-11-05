// Relabelling.cpp
// created by Kuangdai on 6-Jun-2016 
// particle relabelling

#include "Relabelling.h"
#include "Quad.h"
#include "SpectralConstants.h"
#include "Geometric3D.h"
#include "XMath.h"
#include "PreloopFFTW.h"

Relabelling::Relabelling(const Quad *quad):
mMyQuad(quad) {
    mStiff_dZ = RDMatXN::Zero(mMyQuad->getNr(), nPE);
    mStiff_dZdR = RDMatXN::Zero(mMyQuad->getNr(), nPE);
    mStiff_dZdT = RDMatXN::Zero(mMyQuad->getNr(), nPE);
    mStiff_dZdZ = RDMatXN::Zero(mMyQuad->getNr(), nPE);
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            mMass_dZ[ipnt] = RDColX::Zero(mMyQuad->getPointNr(ipol, jpol));
            mMass_dZdR[ipnt] = RDColX::Zero(mMyQuad->getPointNr(ipol, jpol));
            mMass_dZdT[ipnt] = RDColX::Zero(mMyQuad->getPointNr(ipol, jpol));
            mMass_dZdZ[ipnt] = RDColX::Zero(mMyQuad->getPointNr(ipol, jpol));
        }
    }
}

void Relabelling::zeroUndulation() {
    mStiff_dZ.setZero();
    mStiff_dZdR.setZero();
    mStiff_dZdT.setZero();
    mStiff_dZdZ.setZero();
    for (int i = 0; i < nPE; i++) {
        mMass_dZ[i].setZero();
        mMass_dZdR[i].setZero();
        mMass_dZdT[i].setZero();
        mMass_dZdZ[i].setZero();
    }
}

bool Relabelling::isZeroStiff() const {
    if (mStiff_dZ.array().abs().maxCoeff() > tinyDouble) return false;
    return true;
}

bool Relabelling::isZeroMass() const {
    for (int i = 0; i < nPE; i++) 
        if (mMass_dZ[i].array().abs().maxCoeff() > tinyDouble) return false;
    return true;
}

void Relabelling::addUndulation(const Geometric3D &g3D, double srcLat, double srcLon, double srcDep) {
    double rElemCenter = mMyQuad->computeCenterRadius();
    int Nr = mMyQuad->getNr();
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mMyQuad->isAxial());
            const RDMatX3 &rtpS = mMyQuad->computeGeocentricGlobal(srcLat, srcLon, srcDep, xieta, Nr);
            for (int alpha = 0; alpha < Nr; alpha++) {
                double r = rtpS(alpha, 0);
                double t = rtpS(alpha, 1);
                double p = rtpS(alpha, 2);
                mStiff_dZ(alpha, ipnt) += g3D.getDeltaR(r, t, p, rElemCenter); 
            }
        }
    }
}

void Relabelling::finishUndulation() {
    checkHmin();
    formGradientUndulation();
    formMassUndulation();
}

void Relabelling::checkHmin() {
    int Nr = mMyQuad->getNr();
    int maxOrder = (Nr + 1) / 2 - 1;
    RDMatXN original = mStiff_dZ;
    for (int order = 0; order <= maxOrder; order++) {
        // gaussian smooth
        for (int ipnt = 0; ipnt < nPntElem; ipnt++) {
            RDColX data = original.col(ipnt);
            XMath::gaussianSmoothing(data, order, 100, true);
            mStiff_dZ.col(ipnt) = data;
        }
        // compute hmin
        RDColX hmin = mMyQuad->getHminSlices();
        double hmin_all = XMath::trigonResampling(5 * Nr, hmin).minCoeff();
        if (hmin_all >= hmin.minCoeff() * .8) {
            // if (order >= 2) {
            //     std::cout << order << " " << Nr << std::endl;
            //     for (int i = 0; i < Nr; i++) 
            //         std::cout << original.col(0)(i) << " " << mStiff_dZ.col(0)(i) << std::endl;
            //     std::cout << std::endl << std::endl;
            // }
            return;
        }
    }
    throw std::runtime_error("Relabelling::checkHmin || Program should not have reached here.");
}

void Relabelling::formGradientUndulation() {
    int Nr = mMyQuad->getNr();
    int Nu = Nr / 2;
    
    // mask Nyquist first
    if (Nr % 2 == 0) {
        for (int ipnt = 0; ipnt < nPntElem; ipnt++) {
            RDColX data_even = mStiff_dZ.col(ipnt);
            RDColX data_odd = XMath::trigonResampling(Nr - 1, data_even);
            mStiff_dZ.col(ipnt) = XMath::trigonResampling(Nr, data_odd);
        }
    }
    
    // FFT R2C
    vec_CDMatPP deltaR_C(Nu + 1, CDMatPP::Zero());
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            PreloopFFTW::getR2C_RMat(Nr) = mStiff_dZ.col(ipnt);
            PreloopFFTW::computeR2C(Nr);
            CDColX &R2C_C = PreloopFFTW::getR2C_CMat(Nr);
            for (int alpha = 0; alpha <= Nu; alpha++) 
                deltaR_C[alpha](ipol, jpol) = R2C_C(alpha);
        }
    }
    
    // gradient
    const ar3_CDMatPP zero_ar3_CDMatPP = {CDMatPP::Zero(), CDMatPP::Zero(), CDMatPP::Zero()};
    vec_ar3_CDMatPP nablaDeltaR_C(Nu + 1, zero_ar3_CDMatPP);
    mMyQuad->computeGradientScalar(deltaR_C, nablaDeltaR_C, Nu);
    
    // FFT C2R
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            CDColX &C2R_C = PreloopFFTW::getC2R_CMat(Nr);
            // phi = T
            for (int alpha = 0; alpha <= Nu; alpha++) 
                C2R_C(alpha) = nablaDeltaR_C[alpha][1](ipol, jpol);
            PreloopFFTW::computeC2R(Nr);
            mStiff_dZdT.col(ipnt) = PreloopFFTW::getC2R_RMat(Nr);
            // s
            for (int alpha = 0; alpha <= Nu; alpha++) 
                C2R_C(alpha) = nablaDeltaR_C[alpha][0](ipol, jpol);
            PreloopFFTW::computeC2R(Nr);
            RDColX drds = PreloopFFTW::getC2R_RMat(Nr);
            // z
            for (int alpha = 0; alpha <= Nu; alpha++) 
                C2R_C(alpha) = nablaDeltaR_C[alpha][2](ipol, jpol);
            PreloopFFTW::computeC2R(Nr);
            RDColX drdz = PreloopFFTW::getC2R_RMat(Nr);
            
            // rotate s-phi-z to RTZ
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mMyQuad->isAxial());
            double theta = XMath::theta(mMyQuad->mapping(xieta));
            mStiff_dZdZ.col(ipnt) =  drdz * cos(theta) + drds * sin(theta);
            mStiff_dZdR.col(ipnt) = -drdz * sin(theta) + drds * cos(theta);
        }
    }
}

void Relabelling::formMassUndulation() {
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            int nr_mass = mMyQuad->getPointNr(ipol, jpol);
            mMass_dZ[ipnt] = XMath::linearResampling(nr_mass, mStiff_dZ.col(ipnt));
            mMass_dZdR[ipnt] = XMath::linearResampling(nr_mass, mStiff_dZdR.col(ipnt));
            mMass_dZdT[ipnt] = XMath::linearResampling(nr_mass, mStiff_dZdT.col(ipnt));
            mMass_dZdZ[ipnt] = XMath::linearResampling(nr_mass, mStiff_dZdZ.col(ipnt));
        }
    }
}

RDMatXN Relabelling::getStiffJacobian() const {
    int n = mMyQuad->getNr();
    RDMatXN J(n, nPE);
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mMyQuad->isAxial());
            double Z = mMyQuad->mapping(xieta).norm();
            int ipnt = ipol * nPntEdge + jpol;
            if (Z == 0.) {
                const RDColX &J22 = RDColX::Ones(n) + mStiff_dZdZ.col(ipnt);
                J.col(ipnt) = J22.schur(J22).schur(J22);
            } else {
                const RDColX &J00 = RDColX::Ones(n) + mStiff_dZ.col(ipnt) / Z;
                const RDColX &J22 = RDColX::Ones(n) + mStiff_dZdZ.col(ipnt);
                J.col(ipnt) = J00.schur(J00).schur(J22);
            }
        }
    }
    if (J.array().minCoeff() <= 0.) 
        throw std::runtime_error("Relabelling::getStiffJacobian || Negative Jacobian.");
    return J;
}

RDMatXN4 Relabelling::getStiffX() const {
    int n = mMyQuad->getNr();
    RDMatXN4 X(n, nPE * 4);
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mMyQuad->isAxial());
            double Z = mMyQuad->mapping(xieta).norm();
            int ipnt = ipol * nPntEdge + jpol;
            RDColX J0(n), J1(n), J2(n), J3(n);
            if (Z == 0.) {
                J0 = RDColX::Ones(n) + mStiff_dZdZ.col(ipnt);
            } else {
                J0 = RDColX::Ones(n) + mStiff_dZ.col(ipnt) / Z;
            }
            J1 = mStiff_dZdR.col(ipnt);
            J2 = mStiff_dZdT.col(ipnt);
            J3 = RDColX::Ones(n) + mStiff_dZdZ.col(ipnt);
            X.col(nPE * 0 + ipnt) = J0.array().pow(-1.).matrix();
            X.col(nPE * 1 + ipnt) = -J1.schur((J0.schur(J3)).array().pow(-1.).matrix());
            X.col(nPE * 2 + ipnt) = -J2.schur((J0.schur(J3)).array().pow(-1.).matrix());
            X.col(nPE * 3 + ipnt) = J3.array().pow(-1.).matrix();
        }
    }
    return X;
}

arPP_RDColX Relabelling::getMassJacobian() const {
    arPP_RDColX J;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int n = mMyQuad->getPointNr(ipol, jpol);
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mMyQuad->isAxial());
            double Z = mMyQuad->mapping(xieta).norm();
            int ipnt = ipol * nPntEdge + jpol;
            if (Z == 0.) {
                const RDColX &J22 = RDColX::Ones(n) + mMass_dZdZ[ipnt];
                J[ipnt] = J22.schur(J22).schur(J22);
            } else {
                const RDColX &J00 = RDColX::Ones(n) + mMass_dZ[ipnt] / Z;
                const RDColX &J22 = RDColX::Ones(n) + mMass_dZdZ[ipnt];
                J[ipnt] = J00.schur(J00).schur(J22);
            }
            if (J[ipnt].array().minCoeff() <= 0.) 
                throw std::runtime_error("Relabelling::getMassJacobian || Negative Jacobian."); 
        }
    }
    return J;
}

RDMatX3 Relabelling::getSFNormalRTZ(int ipol, int jpol) const {
    int n = mMyQuad->getPointNr(ipol, jpol);
    const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mMyQuad->isAxial());
    double Z = mMyQuad->mapping(xieta).norm();
    int ipnt = ipol * nPntEdge + jpol;
    RDColX J0(n), J1(n), J2(n);
    if (Z == 0.) {
        J0 = RDColX::Ones(n) + mMass_dZdZ[ipnt];
    } else {
        J0 = RDColX::Ones(n) + mMass_dZ[ipnt] / Z;    
    }
    J1 = mMass_dZdR[ipnt];
    J2 = mMass_dZdT[ipnt];
    RDMatX3 normal(n, 3);
    normal.col(0) = -J1.schur(J0);
    normal.col(1) = -J2.schur(J0);
    normal.col(2) = J0.schur(J0);
    return normal;
}

// void Relabelling::dtdR_dtdT(double R, double T, double st, double &dtdR, double &dtdT) {
//     double denominator = sqrt(1. - pow(cos(R) * cos(st) - cos(T) * sin(R) * sin(st), 2.));
//     int iter = 0;
//     while (denominator == 0.) {
//         // singularity
//         R  += (rand() % 201 - 100) * 1e-6;
//         T  += (rand() % 201 - 100) * 1e-6;
//         st += (rand() % 201 - 100) * 1e-6;
//         denominator = sqrt(1. - pow(cos(R) * cos(st) - cos(T) * sin(R) * sin(st), 2.));
//         if (++iter > 100) throw std::runtime_error("Relabelling::dtdR_dtdT || "
//             "Trapped Singularity.");
//     }
//     dtdR = (cos(st) * sin(R) + cos(R) * cos(T) * sin(st)) / denominator;
//     dtdT = -sin(st) * sin(T) / denominator;
// }
// 
// void Relabelling::dpdR_dpdT(double R, double T, double st, double &dpdR, double &dpdT) {
//     double denominator = pow(cos(st) * cos(T) * sin(R) + cos(R) * sin(st), 2.) + pow(sin(R), 2.) * pow(sin(T), 2.);
//     int iter = 0;
//     while (denominator == 0.) {
//         // singularity
//         R  += (rand() % 201 - 100) * 1e-6;
//         T  += (rand() % 201 - 100) * 1e-6;
//         st += (rand() % 201 - 100) * 1e-6;
//         denominator = pow(cos(st) * cos(T) * sin(R) + cos(R) * sin(st), 2.) + pow(sin(R), 2.) * pow(sin(T), 2.);
//         if (++iter > 100) throw std::runtime_error("Relabelling::dpdR_dpdT || Trapped Singularity.");
//     }
//     dpdR = sin(st) * sin(T) / denominator;
//     dpdT = (cos(st) * sin(R) + cos(R) * cos(T) * sin(st)) / denominator;
// }


