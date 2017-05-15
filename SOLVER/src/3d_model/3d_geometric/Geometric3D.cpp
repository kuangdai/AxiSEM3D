// Geometric3D.cpp
// created by Kuangdai on 4-Jun-2016 
// base class of geometric 3D models, such as 
// topography, ellipticity and undulating Moho

#include "Geometric3D.h"

#include "XMPI.h"
#include "XMath.h"
#include "Parameters.h"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

/////////////////////////////// user-defined models here
#include "Ellipticity.h"
#include "Geometric3D_crust1.h"
#include "Geometric3D_Internal.h"
/////////////////////////////// user-defined models here

void Geometric3D::buildInparam(std::vector<Geometric3D *> &models, 
                                const Parameters &par, int verbose) {
    // clear the container
    for (const auto &m: models) {
        delete m;    
    }    
    models.clear();
    
    // check size
    std::size_t nmodels = par.getValue<std::size_t>("MODEL_3D_GEOMETRIC_NUM");
    std::size_t nsize = par.getSize("MODEL_3D_GEOMETRIC_LIST");
    if (nmodels > nsize) {
        throw std::runtime_error("Geometric3D::buildInparam || "
            "Not enough model names provided in MODEL_3D_GEOMETRIC_LIST ||"
            "MODEL_3D_GEOMETRIC_NUM = " + boost::lexical_cast<std::string>(nmodels) + 
            ", but only " + boost::lexical_cast<std::string>(nsize) + " are provided.");
    }
    
    for (std::size_t imodel = 0; imodel < nmodels; imodel++) {
        
        // split model name and parameters
        std::string mstr = par.getValue<std::string>("MODEL_3D_GEOMETRIC_LIST", imodel);
        boost::trim_if(mstr, boost::is_any_of("\t "));
        std::vector<std::string> strs;
        boost::split(strs, mstr, boost::is_any_of("$"), boost::token_compress_on);
        std::string name(strs[0]);
        std::vector<std::string> params(strs.begin() + 1, strs.end());
        
        // create model
        Geometric3D *m;
        if (boost::iequals(name, "crust1")) {
            m = new Geometric3D_crust1();
            
        } else if (boost::iequals(name, "internal")) {
            m = new Geometric3D_Internal();
            
        /////////////////////////////// 
        // user-defined models here
        /////////////////////////////// 
            
        } else {
            throw std::runtime_error("Geometric3D::buildInparam || "
                "Unknown 3D geometric model: " + name + ".");
        }
        
        // initialize
        m->setROuter(XMath::getROuter());
        m->initialize(params);
        models.push_back(m);
        
        // verbose
        if (verbose) {
            XMPI::cout << m->verbose();
        }
    }
    
    // ELLIPTICITY
    std::string emode = par.getValue<std::string>("MODEL_3D_ELLIPTICITY_MODE");
    if (boost::iequals(emode, "full")) {
        Geometric3D *m = new Ellipticity();
        models.push_back(m);
        if (verbose) {
            XMPI::cout << m->verbose();
        }
    }
}

