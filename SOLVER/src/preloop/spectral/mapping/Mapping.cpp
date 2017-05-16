// Mapping.cpp
// created by Kuangdai on 3-May-2016 
// base class of geometric mapping

#include "Mapping.h"

double Mapping::detJacobian(const RDMat24 &nodes, const RDCol2 &xieta, int curvedOuter) const {
    return jacobian(nodes, xieta, curvedOuter).determinant();
}

RDMat22 Mapping::invJacobian(const RDMat24 &nodes, const RDCol2 &xieta, int curvedOuter) const {
    return jacobian(nodes, xieta, curvedOuter).inverse();
}

bool Mapping::invMapping(const RDMat24 &nodes, const RDCol2 &sz, int curvedOuter, RDCol2 &xieta) const {
    xieta = RDCol2::Zero();
    int numiter = 10;
    for (int i = 1; i <= numiter; i++) {
        const RDCol2 &dsz = sz - mapping(nodes, xieta, curvedOuter);
        if (dsz.norm() < 1e-7) {
            return true;
        }
        xieta += invJacobian(nodes, xieta, curvedOuter) * dsz;
    }
    return false;
}

double Mapping::interpolate(const RDRow4 &nodalValues, const RDCol2 &xieta) {
    return interpolateCol(nodalValues, xieta)(0);
}

RDColX Mapping::interpolateCol(const RDMatX4 &nodalValues, const RDCol2 &xieta) {
    RDRow4 shp;
    double xip = 1. + xieta(0);
    double xim = 1. - xieta(0);
    double etap = 1. + xieta(1);
    double etam = 1. - xieta(1);
    shp(0) = xim * etam / 4.;
    shp(1) = xip * etam / 4.;
    shp(2) = xip * etap / 4.;
    shp(3) = xim * etap / 4.;
    return nodalValues * shp.transpose();
}

int Mapping::period0123(int p) {
    if (p > 3) {
        return period0123(p - 4);
    }
    if (p < 0) {
        return period0123(p + 4);
    }
    return p;
}

std::array<RDMat22, 4> Mapping::sOrthogQ2 = {
    (RDMat22(2, 2) << -1., 0., 0., -1.).finished(), 
    (RDMat22(2, 2) << 0., -1., 1., 0.).finished(), 
    (RDMat22(2, 2) << 1., 0., 0., 1.).finished(), 
    (RDMat22(2, 2) << 0., 1., -1., 0.).finished()
};
