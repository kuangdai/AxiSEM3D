// AttBuilder.h
// created by Kuangdai on 14-May-2016 
// compute attenuation parameters 

#pragma once

#include "eigenp.h"

class Parameters;
class Attenuation1D;
class Attenuation3D;
class Quad;

// attenuation parameters from exodus
struct AttParameters {
    AttParameters(int nsls, double fmin, double fmax, double fref, 
        const RDColX &w, const RDColX &y):
        mNSLS(nsls), mFmin(fmin), mFmax(fmax), mFref(fref), mW(w), mY(y) {};
    int mNSLS;
    double mFmin;
    double mFmax;
    double mFref;
    RDColX mW;
    RDColX mY;
};

class AttBuilder {
public:
    AttBuilder(bool cg4, int nsls, double fmin, double fmax, double fref, double deltat): 
        mUseCG4(cg4), mNSLS(nsls), mFmin(fmin), mFmax(fmax), mFref(fref), mDeltaT(deltat) {};
    
    virtual ~AttBuilder() {};
    
    Attenuation1D *createAttenuation1D(const RDMatXN &QKp, const RDMatXN &QMu, 
        RDMatXN &kp, RDMatXN &mu, const Quad *quad) const;
    Attenuation3D *createAttenuation3D(const RDMatXN &QKp, const RDMatXN &QMu, 
        RDMatXN &kp, RDMatXN &mu, const Quad *quad) const;
    std::string verbose() const;
    
    static void buildInparam(AttBuilder *&attBuild, const Parameters &par, 
        const AttParameters *attPar, double dt, int verbose);

protected:    
    virtual void computeAttFactors(const RDMatXN &QKp, const RDMatXN &QMu,
        RDColX &alpha, RDColX &beta, RDColX &gamma,
        RDMatXN &dKpFact, RDMatXN &kpFactAtt, RDMatXN &kpFactNoAtt, 
        RDMatXN &dMuFact, RDMatXN &muFactAtt, RDMatXN &muFactNoAtt) const = 0;
    
    // only the SPECFEM one returns true;    
    virtual bool legacy() const {return false;};
    virtual bool doKappa() const = 0;
    
protected:
    bool mUseCG4;
    int mNSLS;
    double mFmin;
    double mFmax;
    double mFref;
    double mDeltaT;
};
