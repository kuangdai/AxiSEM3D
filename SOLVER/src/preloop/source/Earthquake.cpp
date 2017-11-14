// Earthquake.cpp
// created by Kuangdai on 8-May-2016 
// axial earthquake source

#include "Earthquake.h"
#include "Quad.h"
#include "SpectralConstants.h"
#include <sstream>

#include "Relabelling.h"

Earthquake::Earthquake(double depth, double lat, double lon,
    double Mrr, double Mtt, double Mpp, 
    double Mrt, double Mrp, double Mtp): Source(depth, lat, lon),
    mMxx(Mtt), mMyy(Mpp), mMzz(Mrr), 
    mMxy(Mtp), mMxz(Mrt), mMyz(Mrp) {
    // nothing
}

void Earthquake::computeSourceFourier(const Quad &myQuad, const RDColP &interpFactZ,
    arPP_CMatX3 &fouriers) const {
    // set zero
    for (int ipnt = 0; ipnt < nPntElem; ipnt++) {
        fouriers[ipnt] = CMatX3::Zero(3, 3);
    }
    // Jacobian on axis
    std::array<RDMat22, nPntEdge> axJ;
    int ipol_src = 0;
    for (int jpol_src = 0; jpol_src <= nPol; jpol_src++) {
        const RDCol2 &xieta = SpectralConstants::getXiEta(ipol_src, jpol_src, true);
        axJ[jpol_src] = myQuad.jacobian(xieta);
        axJ[jpol_src] /= axJ[jpol_src].determinant();
    }
    // particle relabelling
    RDColP VX0, VX1, VX2, VX3;
    if (myQuad.hasRelabelling()) {
        const RDMatXN4 &X = myQuad.getRelabelling().getStiffX();
        VX0 = X.block(0, nPE * 0, 1, nPntEdge).transpose();
        VX3 = X.block(0, nPE * 3, 1, nPntEdge).transpose();
        // VX1 and VX2 are phi-independent. The following lines assume phi = 0
        VX1 = X.block(0, nPE * 1, 1, nPntEdge).transpose();
        VX2 = X.block(0, nPE * 2, 1, nPntEdge).transpose();
    } 
    // compute source pointwise
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            // spatial delta function 
            RDMatPP w = RDMatPP::Zero();
            w(ipol, jpol) = 1.;
            const RDMatPP &GU = SpectralConstants::getG_GLJ().transpose() * w;
            const RDMatPP &UG = w * SpectralConstants::getG_GLL();
            for (int jpol_src = 0; jpol_src <= nPol; jpol_src++) {
                double fact = interpFactZ(jpol_src);
                const RDMat22 &J = axJ[jpol_src];
                double dwds = J(1, 1) * GU(ipol_src, jpol_src) - J(1, 0) * UG(ipol_src, jpol_src);
                double dwdz = J(0, 0) * UG(ipol_src, jpol_src) - J(0, 1) * GU(ipol_src, jpol_src);
                if (myQuad.hasRelabelling()) {
                    double X0 = VX0(jpol_src);
                    double X1 = VX1(jpol_src);
                    double X2 = VX2(jpol_src);
                    double X3 = VX3(jpol_src);
                    // monopole
                    fouriers[ipnt](0, 0) += fact * dwds * (mMxx + mMyy) * X0 / (2. * pi);
                    fouriers[ipnt](0, 2) += fact * dwdz * (mMxz * X1 + mMyz * X2 + mMzz * X3) / (2. * pi); 
                    // dipole
                    fouriers[ipnt](1, 0) += fact * dwdz * ((mMxx - iid * mMxy) * X1 + (mMxy - iid * mMyy) * X2 + (mMxz - iid * mMyz) * X3) / (4. * pi);
                    fouriers[ipnt](1, 1) += fact * dwdz * ((mMxx - iid * mMxy) * X1 + (mMxy - iid * mMyy) * X2 + (mMxz - iid * mMyz) * X3) / (4. * pi) * iid;
                    fouriers[ipnt](1, 2) += fact * dwds * (mMxz - iid * mMyz) * X0 / (2. * pi);
                    // quadrupole
                    fouriers[ipnt](2, 0) += fact * dwds * ((mMxx - mMyy) / 2. - iid * mMxy) * X0 / (2. * pi);
                    fouriers[ipnt](2, 1) += fact * dwds * ((mMxx - mMyy) / 2. - iid * mMxy) * X0 / (2. * pi) * iid;
                } else {
                    // monopole
                    fouriers[ipnt](0, 0) += fact * dwds * (mMxx + mMyy) / (2. * pi);
                    fouriers[ipnt](0, 2) += fact * dwdz * mMzz / (2. * pi); 
                    // dipole
                    fouriers[ipnt](1, 0) += fact * dwdz * (mMxz - iid * mMyz) / (4. * pi);
                    fouriers[ipnt](1, 1) += fact * dwdz * (mMxz - iid * mMyz) / (4. * pi) * iid;
                    fouriers[ipnt](1, 2) += fact * dwds * (mMxz - iid * mMyz) / (2. * pi);
                    // quadrupole
                    fouriers[ipnt](2, 0) += fact * dwds * ((mMxx - mMyy) / 2. - iid * mMxy) / (2. * pi);
                    fouriers[ipnt](2, 1) += fact * dwds * ((mMxx - mMyy) / 2. - iid * mMxy) / (2. * pi) * iid;
                }
            }
        }
    }
}

std::string Earthquake::verbose() const {
    std::stringstream ss;
    ss << "\n========================== Source ==========================" << std::endl;
    ss << "  Type         =   " << "Earthquake" << std::endl;
    ss << "  Latitude     =   " << mLatitude << std::endl;
    ss << "  Longitude    =   " << mLongitude << std::endl;
    ss << "  Depth (km)   =   " << mDepth / 1e3 << std::endl;
    ss << "  Moment (N.m) " << std::endl;
    ss << "         Mrr   =   " << (mMzz >= 0. ? " " : "") << mMzz << std::endl;
    ss << "         Mtt   =   " << (mMxx >= 0. ? " " : "") << mMxx << std::endl;
    ss << "         Mpp   =   " << (mMyy >= 0. ? " " : "") << mMyy << std::endl;
    ss << "         Mrt   =   " << (mMxz >= 0. ? " " : "") << mMxz << std::endl;
    ss << "         Mrp   =   " << (mMyz >= 0. ? " " : "") << mMyz << std::endl;
    ss << "         Mtp   =   " << (mMxy >= 0. ? " " : "") << mMxy << std::endl;
    ss << "========================== Source ==========================\n" << std::endl;
    return ss.str();
}



