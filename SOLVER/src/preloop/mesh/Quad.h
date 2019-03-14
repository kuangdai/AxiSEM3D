// Quad.h
// created by Kuangdai on 5-May-2016 
// general Quad element 
// created from Exodus model, elements of AxiSEM3D mesh

#pragma once

#include <vector>
#include "eigenp.h"
#include "Mapping.h"

class ExodusModel;
class NrField;

class Volumetric3D;
class Geometric3D;
class OceanLoad3D;

class GLLPoint;

class Gradient;
class Relabelling;
class Material;
class AttBuilder;

class Domain;

class Quad {
    
public:
    Quad(const ExodusModel &exModel, int quadTag, const NrField &nrf);
    ~Quad();
    
    // 3D models
    void addVolumetric3D(const std::vector<Volumetric3D *> &m3D, 
        double srcLat, double srcLon, double srcDep, double phi2D);
    void addGeometric3D(const std::vector<Geometric3D *> &g3D, 
        double srcLat, double srcLon, double srcDep, double phi2D);
    void setOceanLoad3D(const OceanLoad3D &o3D, 
        double srcLat, double srcLon, double srcDep, double phi2D);
    
    // has relabelling or not
    bool hasRelabelling() const;
    
    // get new deltaT
    double getDeltaT() const;
    
    // setup gll points 
    void setupGLLPoints(std::vector<GLLPoint *> &gllPoints, 
        const IMatPP &myPointTags, double distTol);
    
    // create elements and push to domain
    int release(Domain &domain, const IMatPP &myPointTags, 
        const AttBuilder *attBuild) const;
    int releaseSolid(Domain &domain, const IMatPP &myPointTags, 
        const AttBuilder *attBuild) const;
    int releaseFluid(Domain &domain, const IMatPP &myPointTags) const;
    
    // mapping interfaces
    RDCol2 mapping(const RDCol2 &xieta) const;
    RDMat22 jacobian(const RDCol2 &xieta) const;
    double detJacobian(const RDCol2 &xieta) const;
    bool invMapping(const RDCol2 &sz, RDCol2 &xieta) const;
    
    // for cg4 only
    RDRow4 computeWeightsCG4() const;
    
    // compute geographic coordinates
    RDMatX3 computeGeocentricGlobal(double srcLat, double srcLon, double srcDep,
        const RDCol2 &xieta, int npnt, double phi2D) const;
    double computeCenterRadius() const;
    
    // get properties
    int getQuadTag() const {return mQuadTag;};
    int getElementTag() const {return mElementTag;};
    void setElementTag(int eleTag) {mElementTag = eleTag;};
    bool isFluid() const {return mIsFluid;};
    bool isAxial() const {return mIsAxial;};
    bool onSFBoundary() const {return mOnSFBoundary;};
    bool onSurface() const {return mOnSurface;};
    int getSurfSide() const {return mSurfaceSide;};
    int getNr() const {return mNr;}
    int getNu() const {return mNr / 2;}
    const int &getPointNr(int ipol, int jpol) const {return mPointNr(ipol, jpol);};
    Mapping::MappingTypes getMappingType() const {return mMapping->getType();};
    const RDRowN &getIntegralFactor() const {return mIntegralFactor;};
    const Relabelling &getRelabelling() const {return *mRelabelling;};
    const RDMat24 &getNodalCoords() const {return mNodalCoords;};
    
    // compute gradient of a scalar field
    void computeGradientScalar(const vec_CDMatPP &u, vec_ar3_CDMatPP &u_i) const;
    
    // get hmin on slices
    RDColX getHminSlices() const;
    
    // get field variables for plots
    RDRowN getUndulationOnSlice(double phi) const;
    RDRowN getMaterialOnSlice(const std::string &parName, int refType, double phi) const;
    
    // get spatial range
    void getSpatialRange(double &s_max, double &s_min, double &z_max, double &z_min) const;
    bool nearMe(double s, double z) const;
        
protected:
        
    // create graident operator
    Gradient *createGraident() const;
    
    // form integral factor    
    void formIntegralFactor();
    
    // form mNr and mPointNr
    void formNrField(const NrField &nrf, double distTol);
    
    // get courant number
    double getCourant() const;
    
    // compute normal
    RDMatX3 computeNormal(int side, int ipol, int jpol) const;
    
    ////////////////////////////////////////////////// Nodal level
    // quad tag in Exodus
    int mQuadTag;
    
    // element tag in Domain
    int mElementTag;
    
    // nodal coordinates and tags
    RDMat24 mNodalCoords;
    IRow4 mGlobalNodeTags;
    
    // geometric mapping
    Mapping *mMapping;
    int mCurvedOuter;
    
    // solid or fluid
    bool mIsFluid;
    
    // axial, sf, surface boundaries
    bool mIsAxial, mOnSFBoundary, mOnSurface;
    int mAxialSide, mSFSide, mSurfaceSide;
    
    ////////////////////////////////////////////////// GLL level
    // dt of Quad, nPol = 1
    double mDeltaTRef;
    
    // Nr field
    int mNr;
    IMatPP mPointNr;
    
    // data needed to generate Nr field
    RDRow4 mNodalAveGLLSpacing;
    IRow4 mNearAxisNodes;
    
    // integral factor, just to avoid repeated computation
    RDRowN mIntegralFactor;
    
    // material
    Material *mMaterial; 
    
    // particle relabelling
    Relabelling *mRelabelling = 0;
    
    // Ocean depth
    arPP_RDColX mOceanDepth;
};


