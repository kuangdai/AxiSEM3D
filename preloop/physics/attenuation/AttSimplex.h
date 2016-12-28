// AttSimplex.h
// created by Kuangdai on 14-May-2016 
// SPECFEM method, used for benchmark

#pragma once

#include "AttBuilder.h"

class AttSimplex: public AttBuilder {
public:
    AttSimplex(bool cg4, int nsls, double fmin, double fmax, double fref, double deltat);
    
    void computeFactors(double QMu, double QKappa,
        RDColX &alpha, RDColX &beta, RDColX &gamma,
        double &dKappaFact, double &dMuFact, 
        double &kappaFactAtt, double &muFactAtt,
        double &kappaFactNoAtt, double &muFactNoAtt, bool &doKappa) const;
    
    bool legacy() const {return true;};
    bool doKappa() const {return false;};
        
private:
    
    // intermediate 
    RDColX mTau_s, mFreqs; 
    double mW_central;
    
    const int mNFreq = 100;
};
