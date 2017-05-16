// SemiSphericalMapping.cpp
// created by Kuangdai on 3-May-2016 
// semi-spherical mapping


#include "SemiSphericalMapping.h"
#include "XMath.h"

RDCol2 SemiSphericalMapping::mapping(const RDMat24 &nodes, const RDCol2 &xieta, int curvedOuter) const {
    // rotate system such that curvedOuter = 2
    const RDMat22 &Q2 = sOrthogQ2[curvedOuter];
    const RDMat24 &nodes2 = Q2 * nodes;
    const RDCol2 &xieta2 = Q2 * xieta;
    // get r and theta
    RDMat24 rtheta2;
    rtheta2.row(0).array() = (nodes2.row(0).array().square() + nodes2.row(1).array().square()).sqrt();
    for (int i = 0; i < 4; i++) {
        rtheta2(1, i) = atan2(nodes2(0, i), nodes2(1, i));
    }
    // copy local variables
    double s0 = nodes2(0, Mapping::period0123(curvedOuter - 2));
    double z0 = nodes2(1, Mapping::period0123(curvedOuter - 2));
    double s1 = nodes2(0, Mapping::period0123(curvedOuter - 1)); 
    double z1 = nodes2(1, Mapping::period0123(curvedOuter - 1)); 
    double t2 = rtheta2(1, Mapping::period0123(curvedOuter - 0));
    double r3 = rtheta2(0, Mapping::period0123(curvedOuter + 1));
    double t3 = rtheta2(1, Mapping::period0123(curvedOuter + 1));   
    double xi = xieta2(0);
    double eta = xieta2(1);
    XMath::makeClose(t2, t3);
    // compute in new system
    RDCol2 sz2;
    sz2(0) = (1. + eta) * r3 / 2. * sin(((1. - xi) * t3 + (1. + xi) * t2) / 2.) 
           + (1. - eta) / 2. * (((1. - xi) * s0 + (1. + xi) * s1) / 2.);
    sz2(1) = (1. + eta) * r3 / 2. * cos(((1. - xi) * t3 + (1. + xi) * t2) / 2.) 
           + (1. - eta) / 2. * (((1. - xi) * z0 + (1. + xi) * z1) / 2.); 
    // rotate back        
    return Q2.transpose() * sz2;
}

RDMat22 SemiSphericalMapping::jacobian(const RDMat24 &nodes, const RDCol2 &xieta, int curvedOuter) const {
    // rotate system such that curvedOuter = 2
    const RDMat22 &Q2 = sOrthogQ2[curvedOuter];
    const RDMat24 &nodes2 = Q2 * nodes;
    const RDCol2 &xieta2 = Q2 * xieta;
    // get r and theta
    RDMat24 rtheta2;
    rtheta2.row(0).array() = (nodes2.row(0).array().square() + nodes2.row(1).array().square()).sqrt();
    for (int i = 0; i < 4; i++) {
        rtheta2(1, i) = atan2(nodes2(0, i), nodes2(1, i));
    }
    // copy local variables
    double s0 = nodes2(0, Mapping::period0123(curvedOuter - 2));
    double z0 = nodes2(1, Mapping::period0123(curvedOuter - 2));
    double s1 = nodes2(0, Mapping::period0123(curvedOuter - 1)); 
    double z1 = nodes2(1, Mapping::period0123(curvedOuter - 1)); 
    double t2 = rtheta2(1, Mapping::period0123(curvedOuter - 0));
    double r3 = rtheta2(0, Mapping::period0123(curvedOuter + 1));
    double t3 = rtheta2(1, Mapping::period0123(curvedOuter + 1));     
    double xi = xieta2(0);
    double eta = xieta2(1);
    XMath::makeClose(t2, t3);
    // compute in new system
    RDMat22 J2;
    J2(0, 0) = ((1. - eta) * (s1 - s0)) / 4. 
        + ((1. + eta) * r3 * (t2 - t3) * cos(((1. - xi) * t3 + (1. + xi) * t2) / 2.)) / 4.;
    
    J2(0, 1) = -((1. - xi) * s0 + (1. + xi) * s1) / 4. + r3 / 2. * sin(((1. - xi) * t3 + (1. + xi) * t2) / 2.);
    
    J2(1, 0) = ((1. - eta) * (z1 - z0)) / 4. 
        - ((1. + eta) * r3 * (t2 - t3) * sin(((1. - xi) * t3 + (1. + xi) * t2) / 2.)) / 4.;
    
    J2(1, 1) = -((1. - xi) * z0 + (1. + xi) * z1) / 4. + r3 / 2. * cos(((1. - xi) * t3 + (1. + xi) * t2) / 2.);
    // rotate back        
    return Q2.transpose() * J2 * Q2;
}


