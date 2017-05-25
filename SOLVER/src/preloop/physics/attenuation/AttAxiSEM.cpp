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

void AttAxiSEM::computeAttFactors(const RDMatXN &QKp, const RDMatXN &QMu,
    RDColX &alpha, RDColX &beta, RDColX &gamma,
    RDMatXN &dKpFact, RDMatXN &kpFactAtt, RDMatXN &kpFactNoAtt, 
    RDMatXN &dMuFact, RDMatXN &muFactAtt, RDMatXN &muFactNoAtt) const {
    // size
    alpha = beta = gamma = RDColX(mNSLS);
    dKpFact = kpFactAtt = kpFactNoAtt = QKp;
    dMuFact = muFactAtt = muFactNoAtt = QMu;
    
    // prepare
    double ysum = mY.sum();
    RDColX ydsum = mY / ysum;
    double w_0 = mFref * (2 * pi); 
    double w_1 = sqrt(mFmin * mFmax) * (2. * pi);
    double fact = 0.;
    for (int i = 0; i < mNSLS; i++) {
        fact += ydsum(i) * mW(i) * mW(i) / (w_1 * w_1 + mW(i) * mW(i));
    }
    
    // alpha, beta, gamma
    RDColX onesSLS = RDColX::Ones(mNSLS);
    alpha.array() = (-mW * mDeltaT).array().exp();
    beta.array() = ((onesSLS - alpha).array() / (mW * mDeltaT).array() - alpha.array()) * ydsum.array();
    gamma.array() = ((alpha - onesSLS).array() / (mW * mDeltaT).array() + onesSLS.array()) * ydsum.array();
    
    // size
    RDMatXN ones = RDMatXN::Ones(QKp.rows(), QKp.cols());
    
    // kappa
    if (mDoKappa) {
        kpFactNoAtt = ones + 2. * log(w_1 / w_0) / pi * QKp.array().pow(-1.0).matrix();
        dKpFact.array() = kpFactNoAtt.array() / (1. / ysum * QKp + (1. - fact) * ones).array();
        kpFactAtt = kpFactNoAtt + dKpFact * fact;
    } else {
        dKpFact.setZero();
        kpFactAtt.fill(1.);
        kpFactNoAtt.fill(1.);
    }
    
    // mu
    muFactNoAtt = ones + 2. * log(w_1 / w_0) / pi * QMu.array().pow(-1.0).matrix();
    dMuFact.array() = muFactNoAtt.array() / (1. / ysum * QMu + (1. - fact) * ones).array();
    muFactAtt = muFactNoAtt + dMuFact * fact;    
}
