// Relabelling.h
// created by Kuangdai on 6-Jun-2016 
// particle relabelling

#pragma once 

#include "eigenp.h"

class Quad;
class Geometric3D;

class Relabelling {
public:
    
    Relabelling(const Quad *quad);
    
    void zeroUndulation();
    bool isZeroStiff() const;
    bool isZeroMass() const;
    
    // add deltaR on mass sampling points
    void addUndulation(const Geometric3D &g3D, double srcLat, double srcLon, double srcDep);
    
    // finishUndulation
    void finishUndulation();
    
    // stiffness
    RDMatXN getStiffJacobian() const;
    RDMatXN4 getStiffX() const;
    
    // mass
    arPP_RDColX getMassJacobian() const;
    
    // solid-fluid
    RDMatX3 getSFNormalRTZ(int ipol, int jpol) const;
    
    // deltaR
    const RDMatXN &getDeltaR() const {return mStiff_dZ;};
    
private:
    // check hmin
    void checkHmin();
    
    // compute gradient of deltaR 
    void formGradientUndulation();
    
    // form deltaR on stiffness sampling points (based on mass)
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

