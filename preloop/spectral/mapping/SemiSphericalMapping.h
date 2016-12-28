// SemiSphericalMapping.h
// created by Kuangdai on 3-May-2016 
// semi-spherical mapping


#pragma once

#include "Mapping.h"

class SemiSphericalMapping: public Mapping {
public:
    
    RDCol2 mapping(const RDMat24 &nodes, const RDCol2 &xieta, int curvedOuter) const;
        
    RDMat22 jacobian(const RDMat24 &nodes, const RDCol2 &xieta, int curvedOuter) const;
    
    MappingTypes getType() const {return MappingTypes::SemiSpherical;};
};
