// AttAxiSEM.cpp
// created by Kuangdai on 26-Aug-2016 
// axisem attenuation model
// Geophys. J. Int. (2014) 199, 1078–1093

#include "AttAxiSEM.h"

AttAxiSEM::AttAxiSEM(bool cg4, int nsls, double fmin, double fmax, double fref,
    const RDColX &w, const RDColX &y, double deltat, bool doKappa):
AttBuilder(cg4, nsls, fmin, fmax, fref, deltat), mW(w), mY(y), mDoKappa(doKappa) {
    // nothing
}

void AttAxiSEM::computeFactors(double QMu, double QKappa,
    RDColX &alpha, RDColX &beta, RDColX &gamma,
    double &dKappaFact, double &dMuFact, 
    double &kappaFactAtt, double &muFactAtt,
    double &kappaFactNoAtt, double &muFactNoAtt, bool &doKappa) const {
        
    // prepare
    double ysum = mY.sum();
    RDColX ydsum = mY / mY.sum();
    double w_0 = mFref * (2 * pi); 
    double w_1 = sqrt(mFmin * mFmax) * (2. * pi);
    double fact = 0.;
    for (int i = 0; i < mNSLS; i++) 
        fact += ydsum(i) * mW(i) * mW(i) / (w_1 * w_1 + mW(i) * mW(i));
    
    // alpha, beta, gamma
    alpha = beta = gamma = RDColX(mNSLS);
    RDColX ones = RDColX::Ones(mNSLS);
    alpha.array() = (-mW * mDeltaT).array().exp();
    beta.array() = ((ones - alpha).array() / (mW * mDeltaT).array() - alpha.array()) * ydsum.array();
    gamma.array() = ((alpha - ones).array() / (mW * mDeltaT).array() + ones.array()) * ydsum.array();
    
    // mu
    muFactNoAtt = 1. + 2. / (pi * QMu) * log(w_1 / w_0);
    dMuFact = muFactNoAtt / (1. / ysum * QMu + 1. - fact);
    muFactAtt = muFactNoAtt + dMuFact * fact;
    
    // kappa
    kappaFactNoAtt = 1. + 2. / (pi * QKappa) * log(w_1 / w_0);
    dKappaFact = kappaFactNoAtt / (1. / ysum * QKappa + 1. - fact);
    kappaFactAtt = kappaFactNoAtt + dKappaFact * fact;
    
    // mask kappa
    doKappa = mDoKappa;
    if (!doKappa) {
        dKappaFact = 0.;
        kappaFactAtt = kappaFactNoAtt = 1.;
    }
}

