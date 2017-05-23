// NrFieldEnhance.cpp
// created by Kuangdai on 2-Feb-2017 
// enhanced nr integer field

#include "NrFieldEnhance.h"

#include "Parameters.h"
#include "XMPI.h"
#include <boost/algorithm/string.hpp>
#include "NrFieldEnhanceCircle.h"

NrFieldEnhance::NrFieldEnhance(int ref, bool decrease):
mReference(ref), mDecrease(decrease) {
    // nothing    
}

void NrFieldEnhance::updateNrAtPoint(const RDCol2 &sz_target, int nr_base, int &nr_cur) const {
    double value = getValue(sz_target);
    int newNr = -1;
    if (mReference == 0) {
        newNr = ceil(value);
    } else if (mReference == 1) {
        newNr = ceil(nr_base * (1. + value));
    } else {
        newNr = ceil(nr_cur * (1. + value));
    }
    if (newNr <= 0) return;
    if (mDecrease || newNr > nr_cur) nr_cur = newNr;
}

void NrFieldEnhance::buildInparam(std::vector<NrFieldEnhance *> &nrf, 
    const Parameters &par, int verbose) {
    
    for (const auto &m: nrf) delete m;     
    nrf.clear();
    
    // first check size
    int nnrfs = par.getValue<int>("NU_ENHANCE_NUM");
    int nsize = par.getSize("NU_ENHANCE_LIST");
    if (nnrfs > nsize) throw std::runtime_error("NrFieldEnhance::buildInparam || "
        "Not enough series provided in NU_ENHANCE_LIST ||"
        "NU_ENHANCE_NUM = " + par.getValue<std::string>("NU_ENHANCE_NUM") + ".");
        
    const std::string source = "NrFieldEnhance::buildInparam";
    
    for (int i = 0; i < nnrfs; i++) {
        // split model name and parameters
        std::string mstr = par.getValue<std::string>("NU_ENHANCE_LIST", i);
        std::vector<std::string> strs;
        boost::trim_if(mstr, boost::is_any_of("\t "));
        boost::split(strs, mstr, boost::is_any_of("$"), boost::token_compress_on);
        std::string name = strs[0];
        std::vector<std::string> params;
        for (int i = 1; i < strs.size(); i++) params.push_back(strs[i]);
        
        // create model
        NrFieldEnhance *m;
        if (boost::iequals(name, "circle")) {
            // circle$5171$30$320$320$180$absolute$false
            double r, theta, diameter, hwhm, value;
            std::string strref;
            int ref;
            bool decrease;
            Parameters::castValue(r, params[0], source); r *= 1e3;
            Parameters::castValue(theta, params[1], source); theta *= degree;
            Parameters::castValue(diameter, params[2], source); diameter *= 1e3;
            Parameters::castValue(hwhm, params[3], source); hwhm *= 1e3;
            Parameters::castValue(value, params[4], source); 
            Parameters::castValue(strref, params[5], source); 
            Parameters::castValue(decrease, params[6], source);
            if (boost::iequals(strref, "abs") || boost::iequals(strref, "absolute")) {
                value = value * 2 + 1;
                m = new NrFieldEnhanceCircle(0, decrease, r, theta, diameter, hwhm, value);
            } else if (boost::iequals(strref, "base") || boost::iequals(strref, "ref_base")) {
                m = new NrFieldEnhanceCircle(1, decrease, r, theta, diameter, hwhm, value);
            } else {
                m = new NrFieldEnhanceCircle(2, decrease, r, theta, diameter, hwhm, value);
            }
            
            /////////////////////////////// 
            // user-defined models here
            /////////////////////////////// 
            
        } else {
            throw std::runtime_error("NrFieldEnhance::buildInparam || "
                "Unknown volumetric model name " + name + ".");
        }
        
        // initialize
        if (verbose == 2) XMPI::cout << m->verbose();
        nrf.push_back(m);
    }
}
