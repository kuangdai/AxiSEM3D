// SphericalMapping.cpp
// created by Kuangdai on 3-May-2016 
// spherical mapping


#include "SphericalMapping.h"
#include "XMath.h"

RDCol2 SphericalMapping::mapping(const RDMat24 &nodes, const RDCol2 &xieta, int curvedOuter) const {
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
    double r0 = rtheta2(0, Mapping::period0123(curvedOuter - 2));
    double r3 = rtheta2(0, Mapping::period0123(curvedOuter + 1));
    double t0 = rtheta2(1, Mapping::period0123(curvedOuter - 2));
    double t1 = rtheta2(1, Mapping::period0123(curvedOuter - 1));
    double t2 = rtheta2(1, Mapping::period0123(curvedOuter - 0));
    double t3 = rtheta2(1, Mapping::period0123(curvedOuter + 1));    
    double xi = xieta2(0);
    double eta = xieta2(1);
    XMath::makeClose(t2, t3);
    XMath::makeClose(t0, t1);
    // compute in new system
    RDCol2 sz2;
    sz2(0) = (1. + eta) * r3 / 2. * sin(((1. - xi) * t3 + (1. + xi) * t2) / 2.) 
           + (1. - eta) * r0 / 2. * sin(((1. - xi) * t0 + (1. + xi) * t1) / 2.);
    sz2(1) = (1. + eta) * r3 / 2. * cos(((1. - xi) * t3 + (1. + xi) * t2) / 2.) 
           + (1. - eta) * r0 / 2. * cos(((1. - xi) * t0 + (1. + xi) * t1) / 2.);  
    // rotate back        
    return Q2.transpose() * sz2;
}

RDMat22 SphericalMapping::jacobian(const RDMat24 &nodes, const RDCol2 &xieta, int curvedOuter) const {
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
    double r0 = rtheta2(0, Mapping::period0123(curvedOuter - 2));
    double r3 = rtheta2(0, Mapping::period0123(curvedOuter + 1));
    double t0 = rtheta2(1, Mapping::period0123(curvedOuter - 2));
    double t1 = rtheta2(1, Mapping::period0123(curvedOuter - 1));
    double t2 = rtheta2(1, Mapping::period0123(curvedOuter - 0));
    double t3 = rtheta2(1, Mapping::period0123(curvedOuter + 1));   
    double xi = xieta2(0);
    double eta = xieta2(1);
    XMath::makeClose(t2, t3);
    XMath::makeClose(t0, t1);
    // compute in new system
    RDMat22 J2;
    J2(0, 0) = (1. + eta) * r3 * (t2 - t3) / 4. * cos(((1. - xi) * t3 + (1. + xi) * t2) / 2.) 
             + (1. - eta) * r0 * (t1 - t0) / 4. * cos(((1. - xi) * t0 + (1. + xi) * t1) / 2.);
    
    J2(0, 1) = .5 * (r3 * sin(((1. - xi) * t3 + (1. + xi) * t2) / 2.) 
                   - r0 * sin(((1. - xi) * t0 + (1. + xi) * t1) / 2.));
    
    J2(1, 0) = - (1. + eta) * r3 * (t2 - t3) / 4. * sin(((1. - xi) * t3 + (1. + xi) * t2) / 2.) 
               - (1. - eta) * r0 * (t1 - t0) / 4. * sin(((1. - xi) * t0 + (1. + xi) * t1) / 2.);
    
    J2(1, 1) = .5 * (r3 * cos(((1. - xi) * t3 + (1. + xi) * t2) / 2.) 
                   - r0 * cos(((1. - xi) * t0 + (1. + xi) * t1) / 2.));
    // rotate back 
    return Q2.transpose() * J2 * Q2;
}


