// OceanLoad3D_const.h
// created by Kuangdai on 8-Oct-2016 
// 1D ocean

#pragma once
#include "OceanLoad3D.h"

class OceanLoad3D_const: public OceanLoad3D {
public:

    void initialize(const std::vector<std::string> &params);
    
    double getOceanDepth(double theta, double phi) const {return mDepth;};
    
    std::string verbose() const;
    
private: 
    // PREM value by default
    double mDepth = 3000.;
};
