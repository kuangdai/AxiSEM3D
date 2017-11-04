// AttBuilder.cpp
// created by Kuangdai on 8-Jun-2016 
// compute attenuation parameters 

#include "AttSimplex.h"
#include "AttAxiSEM.h"
#include "Parameters.h"
#include "XMPI.h"
#include <boost/algorithm/string.hpp>
#include "Attenuation1D_Full.h"
#include "Attenuation1D_CG4.h"
#include "Attenuation3D_Full.h"
#include "Attenuation3D_CG4.h"
#include "Quad.h"
#include "XMath.h"

Attenuation1D *AttBuilder::createAttenuation1D(const RDMatXN &QKp, const RDMatXN &QMu, 
    RDMatXN &kp, RDMatXN &mu, const Quad *quad) const {
    // factors
    RDColX alpha, beta, gamma;
    RDMatXN dKpFact, kpFactAtt, kpFactNoAtt;
    RDMatXN dMuFact, muFactAtt, muFactNoAtt;
    computeAttFactors(QKp, QMu,
        alpha, beta, gamma,
        dKpFact, kpFactAtt, kpFactNoAtt, 
        dMuFact, muFactAtt, muFactNoAtt);
        
    // flatten to structured
    RDMatPP kpStr, muStr;
    RDMatPP dKpFactStr, kpFactAttStr, kpFactNoAttStr;
    RDMatPP dMuFactStr, muFactAttStr, muFactNoAttStr;
    XMath::structuredUseFirstRow(kp, kpStr);
    XMath::structuredUseFirstRow(mu, muStr);
    XMath::structuredUseFirstRow(dKpFact, dKpFactStr);
    XMath::structuredUseFirstRow(dMuFact, dMuFactStr);
    XMath::structuredUseFirstRow(kpFactAtt, kpFactAttStr);
    XMath::structuredUseFirstRow(muFactAtt, muFactAttStr);
    XMath::structuredUseFirstRow(kpFactNoAtt, kpFactNoAttStr);
    XMath::structuredUseFirstRow(muFactNoAtt, muFactNoAttStr);
    
    // create attenuation model and change kp and mu in-place
    Attenuation1D *attPtr = 0;
    if (mUseCG4) {
        const RDRow4 &weights_cg4 = quad->computeWeightsCG4();
        RDRow4 dkp, dmu;
        std::vector<int> ipol = {1, 1, 3, 3};
        std::vector<int> jpol = {1, 3, 1, 3};
        for (int i = 0; i < 4; i++) {
            int ip = ipol[i];
            int jp = jpol[i];
            dkp(i) = weights_cg4(i) * dKpFactStr(ip, jp) * kpStr(ip, jp);
            dmu(i) = weights_cg4(i) * dMuFactStr(ip, jp) * muStr(ip, jp);
        }
        kpStr.array() *= kpFactNoAttStr.array();
        muStr.array() *= muFactNoAttStr.array();
        for (int i = 0; i < 4; i++) {
            int ip = ipol[i];
            int jp = jpol[i];
            kpStr(ip, jp) *= 1. + weights_cg4(i) * (kpFactAttStr(ip, jp) / kpFactNoAttStr(ip, jp) - 1.);
            muStr(ip, jp) *= 1. + weights_cg4(i) * (muFactAttStr(ip, jp) / muFactNoAttStr(ip, jp) - 1.);
        }
        attPtr = new Attenuation1D_CG4(mNSLS, 
            alpha.cast<Real>(), beta.cast<Real>(), gamma.cast<Real>(),
            quad->getNu(), dkp.cast<Real>(), dmu.cast<Real>(), doKappa());
    } else {
        RDMatPP dkpStr = dKpFactStr.schur(kpStr);
        RDMatPP dmuStr = dMuFactStr.schur(muStr);
        kpStr.array() *= kpFactAttStr.array();
        muStr.array() *= muFactAttStr.array();
        attPtr = new Attenuation1D_Full(mNSLS, 
            alpha.cast<Real>(), beta.cast<Real>(), gamma.cast<Real>(),
            quad->getNu(), dkpStr.cast<Real>(), dmuStr.cast<Real>(), doKappa());
    }
    
    // structured to flatten
    XMath::flattenFillWithFirstRow(kpStr, kp);
    XMath::flattenFillWithFirstRow(muStr, mu);
    
    return attPtr;
}

Attenuation3D *AttBuilder::createAttenuation3D(const RDMatXN &QKp, const RDMatXN &QMu, 
    RDMatXN &kp, RDMatXN &mu, const Quad *quad) const {
    // factors
    RDColX alpha, beta, gamma;
    RDMatXN dKpFact, kpFactAtt, kpFactNoAtt;
    RDMatXN dMuFact, muFactAtt, muFactNoAtt;
    computeAttFactors(QKp, QMu,
        alpha, beta, gamma,
        dKpFact, kpFactAtt, kpFactNoAtt, 
        dMuFact, muFactAtt, muFactNoAtt);
    
    // create attenuation model and change kp and mu in-place 
    if (mUseCG4) {
        const RDRow4 &weights_cg4 = quad->computeWeightsCG4();
        int nr = quad->getNr();
        RDMatX4 dkp(nr, 4), dmu(nr, 4);
        std::vector<int> ipol = {1, 1, 3, 3};
        std::vector<int> jpol = {1, 3, 1, 3};
        for (int i = 0; i < 4; i++) {
            int ip = ipol[i];
            int jp = jpol[i];
            dkp.col(i) = weights_cg4(i) * dKpFact.col(nPntEdge * ip + jp).schur(kp.col(nPntEdge * ip + jp));
            dmu.col(i) = weights_cg4(i) * dMuFact.col(nPntEdge * ip + jp).schur(mu.col(nPntEdge * ip + jp));
        }
        kp.array() *= kpFactNoAtt.array();
        mu.array() *= muFactNoAtt.array();
        RDColX ones = RDColX::Ones(nr);
        for (int i = 0; i < 4; i++) {
            int ip = ipol[i];
            int jp = jpol[i];
            kp.col(nPntEdge * ip + jp).array() *= ones.array() + weights_cg4(i) * (kpFactAtt.col(nPntEdge * ip + jp).array() 
                / kpFactNoAtt.col(nPntEdge * ip + jp).array() - ones.array());
            mu.col(nPntEdge * ip + jp).array() *= ones.array() + weights_cg4(i) * (muFactAtt.col(nPntEdge * ip + jp).array() 
                / muFactNoAtt.col(nPntEdge * ip + jp).array() - ones.array());
        }
        return new Attenuation3D_CG4(mNSLS, 
            alpha.cast<Real>(), beta.cast<Real>(), gamma.cast<Real>(),
            dkp.cast<Real>(), dmu.cast<Real>(), doKappa());
    } else {
        RDMatXN dkp = dKpFact.schur(kp);
        RDMatXN dmu = dMuFact.schur(mu);
        kp.array() *= kpFactAtt.array();
        mu.array() *= muFactAtt.array();
        return new Attenuation3D_Full(mNSLS, 
            alpha.cast<Real>(), beta.cast<Real>(), gamma.cast<Real>(),
            dkp.cast<Real>(), dmu.cast<Real>(), doKappa());
    }
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

void AttBuilder::buildInparam(AttBuilder *&attBuild, const Parameters &par, 
    const AttParameters *attPar, double dt, int verbose) {
    if (attBuild) {
        delete attBuild;
    }
    
    // pure elastic
    if (!par.getValue<bool>("ATTENUATION") || attPar == 0) {
        attBuild = 0;
        if (verbose) {
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
    
    if (cg4 && nPol != 4) {
        throw std::runtime_error("AttBuilder::buildInparam || "
            "ATTENUATION_CG4 is turned on but nPol is not 4.");
    }
    
    // create model
    if (useLegacy) {
        if (dokappa) {
            throw std::runtime_error("AttBuilder::buildInparam || "
                "Cannot include Q_Kappa when ATTENUATION_SPECFEM_LEGACY = true.");
        }
        if (attPar->mNSLS != 3) {
            throw std::runtime_error("AttBuilder::buildInparam || "
                "Number of SLSs can only be 3 when ATTENUATION_SPECFEM_LEGACY = true.");
        }
        attBuild = new AttSimplex(cg4, attPar->mNSLS, attPar->mFmin, attPar->mFmax, attPar->mFref, dt);
    } else {
        attBuild = new AttAxiSEM(cg4, attPar->mNSLS, attPar->mFmin, attPar->mFmax, attPar->mFref, 
            attPar->mW, attPar->mY, dt, dokappa);
    }
    
    // verbose 
    if (verbose) {
        XMPI::cout << attBuild->verbose();
    } 
}
