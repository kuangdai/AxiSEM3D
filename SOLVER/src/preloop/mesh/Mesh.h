// Mesh.h
// created by Kuangdai on 5-May-2016 
// AxiSEM3D Mesh

#pragma once

#include <vector>
#include "eigenp.h"

class Parameters;
class ExodusModel;
class NrField;
class Volumetric3D;
class Geometric3D;
class OceanLoad3D;
class AttBuilder;

class Quad;
class GLLPoint;
class Domain;
struct DecomposeOption;
struct MessagingInfo;
struct LearnParameters;
class SlicePlot;

class Mesh {
    friend class SlicePlot;
public:
    ~Mesh();
    
    // step 1: define mesh
    Mesh(const ExodusModel *exModel, const NrField *nrf, 
        double srcLat, double srcLon, double srcDep, const Parameters &par, int verbose);
    void setVolumetric3D(const std::vector<Volumetric3D *> &v3D) {mVolumetric3D = v3D;};
    void setGeometric3D(const std::vector<Geometric3D *> &g3D) {mGeometric3D = g3D;};
    void setOceanLoad3D(const OceanLoad3D *o3D) {mOceanLoad3D = o3D;};
    
    // step 2: build unweighted mesh 
    void buildUnweighted();
    
    // step 3: get dt for attenuation
    double getDeltaT() const;
    void setAttBuilder(const AttBuilder *attBuild);
    
    // step 4: build weighted mesh 
    void buildWeighted();
    
    // step 5: release to domain 
    void release(Domain &domain);    
    
    // optional step: test stiffness and mass
    void test();
    
    // get Quads 
    int getNumQuads() const {return mQuads.size();};
    const Quad *getQuad(int index) const {return mQuads[index];};
    
    // solve r such that
    // deltaR(router) + router - depth = deltaR(r) + r
    double computeRadiusRef(double depth, double lat, double lon) const;
    double computeRPhysical(double r, double theta, double phi) const;
    
    // get spatial ranges
    double sMax() const {return mSMax;};
    double sMin() const {return mSMin;};
    double zMax() const {return mZMax;};
    double zMin() const {return mZMin;};
    
    // get max. Nr to initialize solver
    int getMaxNr() const;
    
private:
    
    // build local
    void buildLocal(const DecomposeOption &option);
    
    // destroy local
    void destroy();
    
    // measure
    void measure(DecomposeOption &measured);
    
private:
    
    /////////////////////// global properties ///////////////////////
    // exodus model
    const ExodusModel *mExModel;
    
    // nr field
    const NrField *mNrField;
    
    // source location 
    double mSrcLat, mSrcLon, mSrcDep;
    
    // Volumetric3D models
    std::vector<Volumetric3D *> mVolumetric3D;
    
    // Geometric3D models
    std::vector<Geometric3D *> mGeometric3D;
    
    // OceanLoad3D model
    const OceanLoad3D *mOceanLoad3D;
    
    // attenuation builder
    const AttBuilder *mAttBuilder;
    
    /////////////////////// local build ///////////////////////
    // Quads
    std::vector<Quad *> mQuads;
    
    // GLL Points
    std::vector<GLLPoint *> mGLLPoints;
    
    // element-to-point mapping 
    std::vector<IMatPP> mLocalElemToGLL;
    
    // message info
    MessagingInfo *mMsgInfo;
    
    // spatial ranges
    double mSMax;
    double mSMin;
    double mZMax;
    double mZMin;
    
    ////////////////// domain decomposition //////////////////
    struct DDParameters {
        DDParameters(const Parameters &par);
        bool mReportMeasure;
        int mProcInterval;
        int mNCutsPerProc;
    } *mDDPar;
    
    ////////////////// wisdom learning //////////////////
    LearnParameters *mLearnPar;
    
    ////////////////// 2D in-plane mode //////////////////
    double mPhi2D;
    
    ////////////////// slice plotss //////////////////
    std::vector<SlicePlot *> mSlicePlots;
};



