// PointForce.cpp
// created by Alex on 8-May-2016
// axial point-force source

#include "PointForce.h"
#include "Quad.h"
#include "eigenc.h"
#include "SpectralConstants.h"
#include <sstream>
#include <fstream>

#include "Relabelling.h"

PointForce::PointForce(double depth, double lat, double lon,
    double f1, double f2, double f3): Source(depth, lat, lon),
    px(f1), py(f2), pz(f3) {
    // nothing
}

void PointForce::computeSourceFourier(const Quad &myQuad, const RDColP &interpFactZ,
    arPP_CMatX3 &fouriers) const {
    // set zero
    for (int ipnt = 0; ipnt < nPntElem; ipnt++) {
        fouriers[ipnt] = CMatX3::Zero(2, 3);
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
    RDColP J_PRT;
    if (myQuad.hasRelabelling()) {
        const RDMatXN &JJ = myQuad.getRelabelling().getStiffJacobian();
        J_PRT = JJ.block(0, 0, 1, nPntEdge).transpose();
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
                    //particle relabelling point force
                    // monopole
                    fouriers[ipnt](0, 2) += w(ipol_src, jpol_src) * fact * pz / (2. * pi) * J_PRT[jpol_src];
                    // dipole
                    fouriers[ipnt](1, 0) += w(ipol_src, jpol_src) * fact * (px - iid * py) / (4. * pi) * J_PRT[jpol_src];
                    fouriers[ipnt](1, 1) += w(ipol_src, jpol_src) * fact * (py + iid * px) / (4. * pi) * J_PRT[jpol_src];
                } else {
                    // monopole
                    fouriers[ipnt](0, 2) += w(ipol_src, jpol_src) * fact * pz / (2. * pi);
                    // dipole
                    fouriers[ipnt](1, 0) += w(ipol_src, jpol_src) * fact * (px - iid * py) / (4. * pi);
                    fouriers[ipnt](1, 1) += w(ipol_src, jpol_src) * fact * (py + iid * px) / (4. * pi);
                }
            }
        }
    }
}

std::string PointForce::verbose() const {
    std::stringstream ss;
    ss << "\n========================== Source ==========================" << std::endl;
    ss << "  Type         =   " << "Point Force" << std::endl;
    ss << "  Latitude     =   " << mLatitude << std::endl;
    ss << "  Longitude    =   " << mLongitude << std::endl;
    ss << "  Depth (km)   =   " << mDepth / 1e3 << std::endl;
    ss << "  Component (N) " << std::endl;
    ss << "         f1   =   " << (px >= 0. ? " " : "") << px << std::endl;
    ss << "         f2   =   " << (py >= 0. ? " " : "") << py << std::endl;
    ss << "         f3   =   " << (pz >= 0. ? " " : "") << pz << std::endl;
    ss << "========================== Source ==========================\n" << std::endl;
    return ss.str();
}
