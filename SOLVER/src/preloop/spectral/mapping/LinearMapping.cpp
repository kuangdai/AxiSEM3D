// LinearMapping.cpp
// created by Kuangdai on 3-May-2016 
// linear mapping


#include "LinearMapping.h"

RDCol2 LinearMapping::mapping(const RDMat24 &nodes, const RDCol2 &xieta, int curvedOuter) const {
    return nodes * shapeFunction(xieta).transpose();
}

RDMat22 LinearMapping::jacobian(const RDMat24 &nodes, const RDCol2 &xieta, int curvedOuter) const {    
    return nodes * shapeDerivatives(xieta).transpose();
}

RDRow4 LinearMapping::shapeFunction(const RDCol2 &xieta) const {
    RDRow4 shp;
    double xip = 1. + xieta(0);
    double xim = 1. - xieta(0);
    double etap = 1. + xieta(1);
    double etam = 1. - xieta(1);
    shp(0) = xim * etam / 4.;
    shp(1) = xip * etam / 4.;
    shp(2) = xip * etap / 4.;
    shp(3) = xim * etap / 4.;
    return shp;
}

RDMat24 LinearMapping::shapeDerivatives(const RDCol2 &xieta) const {
    RDMat24 dshp;
    double xip = 1. + xieta(0);
    double xim = 1. - xieta(0);
    double etap = 1. + xieta(1);
    double etam = 1. - xieta(1);
    dshp(0, 0) = - etam / 4.;
    dshp(0, 1) =   etam / 4.;
    dshp(0, 2) =   etap / 4.;
    dshp(0, 3) = - etap / 4.;
    dshp(1, 0) = - xim / 4.;
    dshp(1, 1) = - xip / 4.;
    dshp(1, 2) =   xip / 4.;
    dshp(1, 3) =   xim / 4.;
    return dshp;
}

