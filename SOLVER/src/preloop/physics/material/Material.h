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
class Attenuation1D;
class Attenuation3D;

#include "eigenp.h"

class Material {
public:
    // construction from Exodus
    Material(const Quad *myQuad, const ExodusModel &exModel);
    
    // add 3D
    void addVolumetric3D(const Volumetric3D &m3D, double srcLat, double srcLon, double srcDep, double phi2D);
        
    // Mass
    arPP_RDColX computeElementalMass() const;
    
    // Acoustic
    Acoustic *createAcoustic() const;
    
    // Elastic
    Elastic *createElastic(const AttBuilder *attBuild) const;
        
    // get v_max to compute dt 
    double getVMaxRef() const;
    RDColX getVMax() const;
    
    // check isotropic
    bool isIsotropic() const;
    
    // check 1D
    bool isStiffness1D() const;
    
    // get properties
    RDMatXN getProperty(const std::string &vname, int refType);
        
private:
    // elastic
    Elastic *createElastic1D(const AttBuilder *attBuild) const;
    Elastic *createElastic3D(const AttBuilder *attBuild) const;
    
    // attenuation
    void makeAttenuation1D(const AttBuilder &attBuild, 
        RDMatPP &kappa, RDMatPP &mu, Attenuation1D *&att) const;
    void makeAttenuation3D(const AttBuilder &attBuild, 
        RDMatXN &kappa, RDMatXN &mu, Attenuation3D *&att) const;
    
    // 1D reference
    RDRow4 mVpv1D, mVph1D;
    RDRow4 mVsv1D, mVsh1D;
    RDRow4 mRho1D;
    RDRow4 mEta;
    RDRow4 mQmu, mQkp;
    
    // 3D
    RDMatXN mVpv3D, mVph3D;
    RDMatXN mVsv3D, mVsh3D;
    RDMatXN mRho3D;
    arPP_RDColX mRhoMass3D;
    
    // Quad host
    const Quad *mMyQuad;
};
