// Mapping.h
// created by Kuangdai on 3-May-2016 
// base class of geometric mapping


#pragma once

#include "eigenp.h"
#include <array>

class Mapping {
public:
    
    virtual ~Mapping() {};
        
    enum MappingTypes {Spherical, Linear, SemiSpherical};
    
    virtual RDCol2 mapping(const RDMat24 &nodes, const RDCol2 &xieta, int curvedOuter) const = 0;
        
    virtual RDMat22 jacobian(const RDMat24 &nodes, const RDCol2 &xieta, int curvedOuter) const = 0;
    
    virtual MappingTypes getType() const = 0;
    
    double detJacobian(const RDMat24 &nodes, const RDCol2 &xieta, int curvedOuter) const;
    
    RDMat22 invJacobian(const RDMat24 &nodes, const RDCol2 &xieta, int curvedOuter) const;
    
    bool invMapping(const RDMat24 &nodes, const RDCol2 &sz, int curvedOuter, RDCol2 &xieta) const;
    
    static double interpolate(const RDRow4 &nodalValues, const RDCol2 &xieta);
    
    static RDColX interpolateCol(const RDMatX4 &nodalValues, const RDCol2 &xieta);
    
    static int period0123(int p);

protected:    
    static std::array<RDMat22, 4> sOrthogQ2;
    
};
