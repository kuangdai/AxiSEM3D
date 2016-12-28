// LinearMapping.h
// created by Kuangdai on 3-May-2016 
// linear mapping


#pragma once

#include "Mapping.h"

class LinearMapping: public Mapping {
public:
    
    RDCol2 mapping(const RDMat24 &nodes, const RDCol2 &xieta, int curvedOuter) const;
        
    RDMat22 jacobian(const RDMat24 &nodes, const RDCol2 &xieta, int curvedOuter) const;
    
    MappingTypes getType() const {return MappingTypes::Linear;};
    
private:
    
    RDRow4 shapeFunction(const RDCol2 &xieta) const;
    
    RDMat24 shapeDerivatives(const RDCol2 &xieta) const;
};
