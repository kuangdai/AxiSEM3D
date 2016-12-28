// AttSimplex.cpp
// created by Kuangdai on 14-May-2016 
// SPECFEM method, used for benchmark

#include "AttSimplex.h"

extern "C" {
    void __simplex_MOD_simplex_fminsearch(double *tau_e, const double *tau_s, const int *nsls, 
        const double *f, const int *nf, const double *Qval, int *itercount, double *tolf, int *err);
};

AttSimplex::AttSimplex(bool cg4, int nsls, double fmin, double fmax, double fref, double deltat):
AttBuilder(cg4, nsls, fmin, fmax, fref, deltat) {
    double log_min = std::log10(mFmin);
    double log_max = std::log10(mFmax);
    double log_sum = log_min + log_max;
    double log_dif = log_max - log_min;  
    mTau_s =  RDColX::Zero(mNSLS);
    for (int i = 0; i < mNSLS; i++) {
        double log_i = log_min + i * log_dif / (mNSLS - 1);
        mTau_s(i) = 1.0 / (2. * pi * pow(10., log_i));
    }
    mFreqs = RDColX::Zero(mNFreq);
    for (int i = 0; i < mNFreq; i++) 
        mFreqs(i) = log_min + i * log_dif / (mNFreq - 1);
    mW_central = 2. * pi * pow(10., log_sum / 2.);
}

void AttSimplex::computeFactors(double QMu, double QKappa,
    RDColX &alpha, RDColX &beta, RDColX &gamma,
    double &dKappaFact, double &dMuFact, 
    double &kappaFactAtt, double &muFactAtt,
    double &kappaFactNoAtt, double &muFactNoAtt, bool &doKappa) const {
    
    // no kappa attenuation
    doKappa = false;
    dKappaFact = 0.;
    kappaFactAtt = kappaFactNoAtt = 1.;
    
    // minimization of frequency misfit
    RDColX tau_e = mTau_s * (1. + 2. / QMu);
    int niter, ierr; double tol;
    __simplex_MOD_simplex_fminsearch(tau_e.data(), mTau_s.data(), &mNSLS, 
        mFreqs.data(), &mNFreq, &QMu, &niter, &tol, &ierr);
    if (ierr > 0) throw std::runtime_error("AttSimplex::computeFactors || "
        "Convergence failed in Fortran subroutine simplex_fminsearch");
    
    // alpha beta gamma
    RDColX betax = RDColX::Ones(mNSLS) - (tau_e.array() / mTau_s.array()).matrix();
    RDColX tauinv = -mTau_s.array().pow(-1).matrix();
    RDColX factor_common = betax.schur(tauinv);
    alpha = (RDColX::Ones(mNSLS).array() + mDeltaT * tauinv.array() 
            + pow(mDeltaT, 2) * tauinv.array().pow(2) / 2. 
            + pow(mDeltaT, 3) * tauinv.array().pow(3) / 6. 
            + pow(mDeltaT, 4) * tauinv.array().pow(4) / 24.).matrix();
    beta = (RDColX::Ones(mNSLS).array() * mDeltaT / 2. 
            + pow(mDeltaT, 2) * tauinv.array() / 3. 
            + pow(mDeltaT, 3) * tauinv.array().pow(2) / 8. 
            + pow(mDeltaT, 4) * tauinv.array().pow(3) / 24.).matrix().schur(factor_common);
    gamma = (RDColX::Ones(mNSLS).array() * mDeltaT / 2. 
            + pow(mDeltaT, 2) * tauinv.array() / 6. 
            + pow(mDeltaT, 3) * tauinv.array().pow(2) / 24.).matrix().schur(factor_common);            
    
    // factors
    double factor_scale_mu0 = 1. + 2. * log(mW_central / (2. * pi)) / (pi * QMu);
    double one_minus_sum_beta = 1. - betax.sum();
    double a_val = 1.;
    double b_val = 0.;
    for (int i = 0; i < mNSLS; i++) {
        double temp = 1. + mW_central * mW_central * tau_e(i) * tau_e(i);
        a_val -= mW_central * mW_central * tau_e(i) * (tau_e(i) - mTau_s(i)) / temp;
        b_val += mW_central * (tau_e(i) - mTau_s(i)) / temp;
    }
    double big_omega = a_val * (sqrt(1. + b_val * b_val / (a_val * a_val)) - 1.);
    double factor_scale_mu = b_val * b_val / (2. * big_omega);
    dMuFact = factor_scale_mu0 * factor_scale_mu;
    muFactAtt = factor_scale_mu0 * factor_scale_mu * one_minus_sum_beta;
    muFactNoAtt = factor_scale_mu0;
}



