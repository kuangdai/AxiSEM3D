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
    Elastic *createElasticAniso(bool elem1D, const AttBuilder *attBuild) const;
        
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
        
    //////////// anisotropy /////////////
private:    
    void initAniso();
    void rotateAniso(double srcLat, double srcLon, double srcDep);
    static RDMatXX bondTransformation(RDMatXX inCijkl, double alpha, double beta, double gamma);
    
private:
    void prepare3D();
    bool _3Dprepared() const;
    
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
    
    // full anisotropy in 1D and 3D
    RDRow4 mC11_1D, mC12_1D, mC13_1D, mC14_1D, mC15_1D, mC16_1D;
    RDRow4 mC22_1D, mC23_1D, mC24_1D, mC25_1D, mC26_1D;
    RDRow4 mC33_1D, mC34_1D, mC35_1D, mC36_1D;
    RDRow4 mC44_1D, mC45_1D, mC46_1D;
    RDRow4 mC55_1D, mC56_1D;
    RDRow4 mC66_1D;
    RDMatXN mC11_3D, mC12_3D, mC13_3D, mC14_3D, mC15_3D, mC16_3D;
    RDMatXN mC22_3D, mC23_3D, mC24_3D, mC25_3D, mC26_3D;
    RDMatXN mC33_3D, mC34_3D, mC35_3D, mC36_3D;
    RDMatXN mC44_3D, mC45_3D, mC46_3D;
    RDMatXN mC55_3D, mC56_3D;
    RDMatXN mC66_3D;
    bool mFullAniso = false;
    
    // 3D material sampled at mass points
    arPP_RDColX mRhoMass3D;
    arPP_RDColX mVpFluid3D;
    
    // Quad host
    const Quad *mMyQuad;
};
