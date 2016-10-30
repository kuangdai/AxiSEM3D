// OceanLoad3D.cpp
// created by Kuangdai on 8-Oct-2016 
// base class of Ocean Load 3D models

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
    if (model) delete model;
    
    std::string mstr = par.getValue<std::string>("MODEL_OCEAN_LOAD");
    if (boost::iequals(mstr, "none")) {
        model = 0;
    } else {
        // split model name and parameters
        std::vector<std::string> strs;
        boost::trim_if(mstr, boost::is_any_of("\t "));
        boost::split(strs, mstr, boost::is_any_of("$"), boost::token_compress_on);
        std::string name = strs[0];
        std::vector<double> params;
        try { 
            for (int i = 1; i < strs.size(); i++) 
                params.push_back(boost::lexical_cast<double>(strs[i]));
        } catch (std::exception) {
            throw std::runtime_error("OceanLoad3D::buildInparam || "
                "Invalid parameter following ocean load model " + name + ".");
        }
        
        if (boost::iequals(name, "constant")) {
            model = new OceanLoad3D_const();
        } else if (boost::iequals(name, "crust1")) {
            model = new OceanLoad3D_crust1();
            
        /////////////////////////////// 
        // user-defined models here
        /////////////////////////////// 
            
        } else {
            throw std::runtime_error("OceanLoad3D::buildInparam || "
                "Unknown ocean load model name " + name + ".");
        }
        
        model->initialize(params);
        if (verbose) XMPI::cout << model->verbose();
    }
}
