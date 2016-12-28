// AttBuilder.h
// created by Kuangdai on 14-May-2016 
// compute attenuation parameters 

#pragma once

#include "eigenp.h"

class Parameters;

struct AttParameters {
    AttParameters(int nsls, double fmin, double fmax, double fref, RDColX w, RDColX y):
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
    
    virtual void computeFactors(double QMu, double QKappa,
        RDColX &alpha, RDColX &beta, RDColX &gamma,
        double &dKappaFact, double &dMuFact, 
        double &kappaFactAtt, double &muFactAtt,
        double &kappaFactNoAtt, double &muFactNoAtt, bool &doKappa) const = 0;
    
    virtual bool legacy() const = 0;
    virtual bool doKappa() const = 0;
            
    int getNSLS() const {return mNSLS;};
    bool useCG4() const {return mUseCG4;};
    
    std::string verbose() const;
        
    static void buildInparam(AttBuilder *&attBuild, const Parameters &par, 
        const AttParameters &attPar, double dt, int verbose);
    
protected:
    bool mUseCG4;
    int mNSLS;
    double mFmin;
    double mFmax;
    double mFref;
    double mDeltaT;
};
