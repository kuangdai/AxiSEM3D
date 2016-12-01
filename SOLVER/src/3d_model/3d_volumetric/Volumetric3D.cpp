// Volumetric3D.cpp
// created by Kuangdai on 16-May-2016 
// base class of Volumetric 3D models, including both mantle and crust
// currently we assume 1D Q-factors (Q_mu and Q_kappa)

#include "Volumetric3D.h"

#include "XMPI.h"
#include "XMath.h"
#include "Parameters.h"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

/////////////////////////////// user-defined models here
#include "Volumetric3D_s20rts.h"
#include "Volumetric3D_s40rts.h"
#include "Volumetric3D_crust1.h"
#include "Volumetric3D_bubble.h"
#include "Volumetric3D_cylinder.h"
/////////////////////////////// user-defined models here

void Volumetric3D::buildInparam(std::vector<Volumetric3D *> &models, 
    const Parameters &par, int verbose) {
    for (const auto &m: models) delete m;     
    models.clear();
    // first check size
    int nmodels = par.getValue<int>("MODEL_3D_VOLUMETRIC_NUM");
    int nsize = par.getSize("MODEL_3D_VOLUMETRIC_LIST");
    if (nmodels > nsize) throw std::runtime_error("Volumetric3D::buildInparam || "
        "Not enough model names provided in MODEL_3D_VOLUMETRIC_LIST ||"
        "MODEL_3D_VOLUMETRIC_NUM = " + par.getValue<std::string>("MODEL_3D_VOLUMETRIC_NUM") + ".");
    
    for (int i = 0; i < nmodels; i++) {
        // split model name and parameters
        std::string mstr = par.getValue<std::string>("MODEL_3D_VOLUMETRIC_LIST", i);
        std::vector<std::string> strs;
        boost::trim_if(mstr, boost::is_any_of("\t "));
        boost::split(strs, mstr, boost::is_any_of("$"), boost::token_compress_on);
        std::string name = strs[0];
        std::vector<std::string> params;
        for (int i = 1; i < strs.size(); i++) params.push_back(strs[i]);
        
        // create model
        Volumetric3D *m;
        if (boost::iequals(name, "s20rts")) {
            m = new Volumetric3D_s20rts();
        } else if (boost::iequals(name, "s40rts")) {
            m = new Volumetric3D_s40rts();
        } else if (boost::iequals(name, "crust1")) {
            m = new Volumetric3D_crust1();
        } else if (boost::iequals(name, "bubble")) {
            m = new Volumetric3D_bubble(); 
        } else if (boost::iequals(name, "cylinder")) {
            m = new Volumetric3D_cylinder();        
            
        /////////////////////////////// 
        // user-defined models here
        /////////////////////////////// 
            
        } else {
            throw std::runtime_error("Volumetric3D::buildInparam || "
                "Unknown volumetric model name " + name + ".");
        }
        
        // initialize
        m->setROuter(XMath::getROuter());
        m->initialize(params);
        if (verbose) XMPI::cout << m->verbose();
        models.push_back(m);
    }
}
