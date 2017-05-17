// OceanLoad3D_const.cpp
// created by Kuangdai on 8-Oct-2016 
// 1D ocean

// Ellipticity.cpp
// created by Kuangdai on 4-Jun-2016 
// ellipticity

#include "OceanLoad3D_const.h"
#include <sstream>
#include "Parameters.h"

void OceanLoad3D_const::initialize(const std::vector<std::string> &params) {
    if (params.size() >= 1) {
        Parameters::castValue(mDepth, params[0], "OceanLoad3D_const::initialize");
        mDepth *= 1e3;
    }
}

std::string OceanLoad3D_const::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D OceanLoad =======================" << std::endl;
    ss << "  Model Name   =   Constant" << std::endl;
    ss << "  Depth / km   =   " << mDepth / 1e3 << std::endl;
    ss << "======================= 3D OceanLoad =======================\n" << std::endl;
    return ss.str();
}

