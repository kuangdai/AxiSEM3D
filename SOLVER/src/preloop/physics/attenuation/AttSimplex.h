// AttSimplex.h
// created by Kuangdai on 14-May-2016 
// SPECFEM method, used for benchmark

#pragma once

#include "AttBuilder.h"

class AttSimplex: public AttBuilder {
public:
    AttSimplex(bool cg4, int nsls, double fmin, double fmax, double fref, double deltat);
    
    void computeAttFactors(const RDMatXN &QKp, const RDMatXN &QMu,
        RDColX &alpha, RDColX &beta, RDColX &gamma,
        RDMatXN &dKpFact, RDMatXN &kpFactAtt, RDMatXN &kpFactNoAtt, 
        RDMatXN &dMuFact, RDMatXN &muFactAtt, RDMatXN &muFactNoAtt) const;
    bool doKappa() const {return false;};
        
private:
    
    // intermediate 
    RDColX mTau_s, mFreqs; 
    double mW_central;
    
    const int mNFreq = 100;
};
