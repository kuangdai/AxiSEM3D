// Material.h
// created by Kuangdai on 17-May-2016 
// 3D seismic material properties

#pragma once

class Quad;
class ExodusModel;
class Volumetric3D;

class Acoustic;
class Elastic;
class AttBuilder;

#include "eigenp.h"
#include <array>

class Material {
public:
    // construction from Exodus
    Material(const Quad *myQuad, const ExodusModel &exModel);
    
    // add 3D
    void addVolumetric3D(const std::vector<Volumetric3D *> &m3D, 
        double srcLat, double srcLon, double srcDep, double phi2D);
        
    // Mass
    arPP_RDColX computeElementalMass() const;
    
    // Acoustic
    Acoustic *createAcoustic(bool elem1D) const;
    
    // Elastic
    Elastic *createElastic(bool elem1D, const AttBuilder *attBuild) const;
        
    // get v_max to compute dt 
    double getVMaxRef() const;
    RDColX getVMax() const;
    
    // check 1D
    bool isFluidPar1D() const;
    bool isSolidPar1D(bool attenuation) const;
    
    // check isotropic
    bool isIsotropic() const;
    
    // get properties
    RDMatXN getProperty(const std::string &vname, int refType);
        
private:
    // 1D reference material
    RDRow4 mVpv1D, mVph1D;
    RDRow4 mVsv1D, mVsh1D;
    RDRow4 mRho1D;
    RDRow4 mEta1D;
    RDRow4 mQkp1D, mQmu1D;
    
    // 3D material
    RDMatXN mVpv3D, mVph3D;
    RDMatXN mVsv3D, mVsh3D;
    RDMatXN mRho3D;
    RDMatXN mEta3D;
    RDMatXN mQkp3D, mQmu3D;
    
    // 3D material sampled at mass points
    arPP_RDColX mRhoMass3D;
    arPP_RDColX mVpFluid3D;
    
    // Quad host
    const Quad *mMyQuad;
};
