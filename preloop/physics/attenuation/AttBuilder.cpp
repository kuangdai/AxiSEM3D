// AttBuilder.cpp
// created by Kuangdai on 8-Jun-2016 
// compute attenuation parameters 

#include "AttSimplex.h"
#include "AttAxiSEM.h"
#include "Parameters.h"
#include "XMPI.h"
#include <boost/algorithm/string.hpp>

void AttBuilder::buildInparam(AttBuilder *&attBuild, const Parameters &par, 
    const AttParameters &attPar, double dt, int verbose) {
    if (attBuild) delete attBuild;
    
    // pure elastic
    if (!par.getValue<bool>("ATTENUATION")) {
        attBuild = 0;
        if (verbose == 2) {
            XMPI::cout << "\n=================== Attenuation Builder ====================" << XMPI::endl;
            XMPI::cout << "  Attenuation is turned off. "  << XMPI::endl;
            XMPI::cout << "=================== Attenuation Builder ====================\n" << XMPI::endl;
        }
        return;
    }
    
    // parameters
    bool useLegacy = par.getValue<bool>("ATTENUATION_SPECFEM_LEGACY");
    bool cg4 = par.getValue<bool>("ATTENUATION_CG4");
    bool dokappa = par.getValue<bool>("ATTENUATION_QKAPPA");
    
    // create model
    if (useLegacy) {
        if (dokappa) throw std::runtime_error("AttBuilder::buildInparam || "
            "Cannot include Q_Kappa when ATTENUATION_SPECFEM_LEGACY = true.");
        attBuild = new AttSimplex(cg4, attPar.mNSLS, attPar.mFmin, attPar.mFmax, attPar.mFref, dt);
    } else {
        attBuild = new AttAxiSEM(cg4, attPar.mNSLS, attPar.mFmin, attPar.mFmax, attPar.mFref, 
            attPar.mW, attPar.mY, dt, dokappa);
    }
    
    // verbose 
    if (verbose == 2) XMPI::cout << attBuild->verbose(); 
}

std::string AttBuilder::verbose() const {
    std::stringstream ss;
    ss << "\n=================== Attenuation Builder ====================" << std::endl;
    ss << "  Coarse Grained    =   " << (mUseCG4 ? "YES" : "NO") << std::endl;
    ss << "  Number of SLS     =   " << mNSLS << std::endl;
    ss << "  Freq. Band (Hz)   =   [" << mFmin << ", " << mFmax << "]" << std::endl;
    ss << "  Ref. Freq. (Hz)   =   " << mFref << std::endl;
    ss << "  Time Step         =   " << mDeltaT << std::endl;
    ss << "  Include QKappa    =   " << (doKappa() ? "YES" : "NO") << std::endl;
    if (legacy()) ss << "  Using SPECFEM Legacy model." << std::endl;
    ss << "=================== Attenuation Builder ====================\n" << std::endl;
    return ss.str();
}
