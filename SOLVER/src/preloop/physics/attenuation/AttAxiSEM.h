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
    
    void computeFactors(double QMu, double QKappa,
        RDColX &alpha, RDColX &beta, RDColX &gamma,
        double &dKappaFact, double &dMuFact, 
        double &kappaFactAtt, double &muFactAtt,
        double &kappaFactNoAtt, double &muFactNoAtt, bool &doKappa) const;
    
    bool legacy() const {return false;};
    bool doKappa() const {return mDoKappa;};
        
private:
    
    RDColX mW, mY;
    bool mDoKappa;
};
