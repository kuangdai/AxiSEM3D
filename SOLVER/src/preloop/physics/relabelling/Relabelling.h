// Relabelling.h
// created by Kuangdai on 6-Jun-2016 
// particle relabelling

#pragma once 

#include "eigenp.h"

class Quad;
class Geometric3D;
class PRT;

class Relabelling {
public:
    
    Relabelling(const Quad *quad);
    
    // add deltaR on mass sampling points
    void addUndulation(const std::vector<Geometric3D *> &g3D, 
        double srcLat, double srcLon, double srcDep, double phi2D);
    
    // stiffness
    RDMatXN getStiffJacobian() const;
    RDMatXN4 getStiffX() const;
    
    // mass
    arPP_RDColX getMassJacobian() const;
    
    // solid-fluid
    RDMatX3 getSFNormalRTZ(int ipol, int jpol) const;
    
    // is undulation zero
    bool isZero() const;
    
    // is undulation 1D
    bool isPar1D() const;
    
    // create PRT pointer
    PRT *createPRT(bool elem1D) const;
    
    // deltaR
    const RDMatXN &getDeltaR() const {return mStiff_dZ;};
    
private:
    // check hmin
    void checkHmin();
    
    // compute gradient of deltaR 
    void formGradientUndulation();
    
    // form deltaR on mass sampling points (based on stiffness)
    void formMassUndulation();
    
    // host Quad 
    const Quad *mMyQuad;
    // undulation derivatives
    // in this class, the variables are defined in RTZ coordinate system
    // R = Radial = theta
    // T = Transverse = phi
    // Z = vertical = r 
    RDMatXN mStiff_dZ;
    RDMatXN mStiff_dZdR;
    RDMatXN mStiff_dZdT;
    RDMatXN mStiff_dZdZ;
    arPP_RDColX mMass_dZ;
    arPP_RDColX mMass_dZdR;
    arPP_RDColX mMass_dZdT;
    arPP_RDColX mMass_dZdZ;
};

