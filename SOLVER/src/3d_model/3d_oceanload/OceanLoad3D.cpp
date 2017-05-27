// OceanLoad3D.cpp
// created by Kuangdai on 8-Oct-2016 
// base class of Ocean Load models

#include "OceanLoad3D.h"

#include "XMPI.h"
#include "Parameters.h"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

/////////////////////////////// user-defined models here
#include "OceanLoad3D_const.h"
#include "OceanLoad3D_crust1.h"
/////////////////////////////// user-defined models here

// change density of ocean water here 
const double OceanLoad3D::mWaterDensity = 1027.;

void OceanLoad3D::buildInparam(OceanLoad3D *&model,
                               const Parameters &par, int verbose) {
    // delete    
    if (model) {
        delete model;
    }
    
    // get keyword
    std::string mstr = par.getValue<std::string>("MODEL_3D_OCEAN_LOAD");
    std::vector<std::string> strs = Parameters::splitString(mstr, "$");
    std::string name(strs[0]);
    std::vector<std::string> params(strs.begin() + 1, strs.end());
    
    // create model
    if (boost::iequals(name, "none")) {
        model = 0;
    } else if (boost::iequals(name, "constant")) {
        model = new OceanLoad3D_const();
    } else if (boost::iequals(name, "crust1")) {
        model = new OceanLoad3D_crust1();
    
    /////////////////////////////// 
    // user-defined models here
    /////////////////////////////// 
        
    } else {
        throw std::runtime_error("OceanLoad3D::buildInparam || "
            "Unknown ocean-load model name " + name + ".");
    }
    
    if (model != 0) {
        // initialize model
        model->initialize(params);
        
        // verbose
        if (verbose) {
            XMPI::cout << model->verbose();
        }
    }
}
