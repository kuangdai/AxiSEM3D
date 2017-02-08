// NrField.cpp
// created by Kuangdai on 13-May-2016 
// base class of nr integer field

#include "NrField.h"
#include "ConstNrField.h"
#include "EmpNrField.h"
#include "WisdomNrField.h"
#include "NrFieldEnhance.h"

#include "Parameters.h"
#include "XMPI.h"
#include <boost/algorithm/string.hpp>

NrField::~NrField() {
    for (const auto &m: mEnhance) delete m;
}

void NrField::buildInparam(NrField *&nrf, const Parameters &par, 
    double router, int verbose) {
    if (nrf) delete nrf;
    std::string type = par.getValue<std::string>("NU_TYPE");
    bool useLucky = par.getValue<bool>("NU_FFTW_LUCKY_NUMBER");
    
    if (boost::iequals(type, "constant")) {
        int nu = par.getValue<int>("NU_CONST");
        nrf = new ConstNrField(useLucky, nu);
    } else if (boost::iequals(type, "empirical")) {
        int nu_ref = par.getValue<int>("NU_EMP_REF");
        int nu_min = par.getValue<int>("NU_EMP_MIN");
        bool scaleS = par.getValue<bool>("NU_EMP_SCALE_AXIS");
        bool scaleT = par.getValue<bool>("NU_EMP_SCALE_THETA");
        bool scaleD = par.getValue<bool>("NU_EMP_SCALE_DEPTH");
        double powS = par.getValue<double>("NU_EMP_POW_AXIS");
        double factPI = par.getValue<double>("NU_EMP_FACTOR_PI");
        double startT = par.getValue<double>("NU_EMP_THETA_START") * degree;
        double powT = par.getValue<double>("NU_EMP_POW_THETA");
        double factD0 = par.getValue<double>("NU_EMP_FACTOR_SURF");
        double startD = par.getValue<double>("NU_EMP_DEPTH_START") * 1e3;
        double endD = par.getValue<double>("NU_EMP_DEPTH_END") * 1e3;
        nrf = new EmpNrField(useLucky, nu_ref, nu_min, scaleS, scaleT, scaleD, 
            router, powS, factPI, startT, powT, factD0, startD, endD);
    } else if (boost::iequals(type, "wisdom")) {
        std::string fname = Parameters::sInputDirectory + "/" + par.getValue<std::string>("NU_WISDOM_REUSE_INPUT");
        double factor = par.getValue<double>("NU_WISDOM_REUSE_FACTOR");
        nrf = new WisdomNrField(useLucky, fname, factor);
    } else {
        throw std::runtime_error("NrField::build || "
            "Invalid parameter, keyword = NU_TYPE.");
    }
    
    if (verbose == 2) XMPI::cout << nrf->verbose();
    
    // enhanced
    NrFieldEnhance::buildInparam(nrf->mEnhance, par, verbose);
}

int NrField::enhancedNr(const RDCol2 &coords, int nr_base) const {
    int nr_cur = nr_base;
    for (const auto &m: mEnhance) m->updateNrAtPoint(coords, nr_base, nr_cur);
    return nr_cur;
}
