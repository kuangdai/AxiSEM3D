// AttAxiSEM.h
// created by Kuangdai on 26-Aug-2016 
// axisem attenuation model
// Geophys. J. Int. (2014) 199, 1078–1093

#pragma once

#include "AttBuilder.h"

class AttAxiSEM: public AttBuilder {
public:
    AttAxiSEM(bool cg4, int nsls, double fmin, double fmax, double fref,
        const RDColX &w, const RDColX &y, double deltat, bool doKappa);
    
    void computeAttFactors(const RDMatXN &QKp, const RDMatXN &QMu,
        RDColX &alpha, RDColX &beta, RDColX &gamma,
        RDMatXN &dKpFact, RDMatXN &kpFactAtt, RDMatXN &kpFactNoAtt, 
        RDMatXN &dMuFact, RDMatXN &muFactAtt, RDMatXN &muFactNoAtt) const;
    bool doKappa() const {return mDoKappa;};
        
private:
    
    RDColX mW, mY;
    bool mDoKappa;
};
