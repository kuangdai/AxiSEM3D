// PointwiseRecorder.h
// created by Kuangdai on 1-Jun-2017 
// recorder for point-wise receivers

#pragma once

#include "eigenc.h"
#include "eigenp.h"
class Element;
class PointwiseIO;

// receiver info
struct PointwiseInfo {
    PointwiseInfo(const std::string &name, const std::string &network, 
        double phi, const RDMatPP &weights, const Element *ele,
        double theta, double baz, 
        double lat, double lon, double dep, bool dumpStrain, bool dumpCurl):
        mName(name), mNetwork(network), 
        mPhi(phi), mWeights(weights.cast<Real>()), mElement(ele),
        mTheta(theta), mBAz(baz),
        mLat(lat), mLon(lon), mDep(dep), mDumpStrain(dumpStrain), mDumpCurl(dumpCurl) {};
    
    //// name
    std::string mName;
    std::string mNetwork;
    
    //// to compute disp from Element
    double mPhi;
    RMatPP mWeights;
    const Element *mElement;
    
    //// to transform disp to RTZ or ENZ
    double mTheta;
    double mBAz;
    
    // for station XML
    double mLat, mLon, mDep;
    
    // dump strain
    bool mDumpStrain;
    
    // dump curl
    bool mDumpCurl;
};

class PointwiseRecorder {
public:
    PointwiseRecorder(int totalRecordSteps, int recordInterval, 
        int bufferSize, const std::string &components, 
        double srcLat, double srcLon, double srcDep);
    ~PointwiseRecorder();
    
    // add a receiver
    void addReceiver(const std::string &name, const std::string &network,
        double phi, const RDMatPP &weights, const Element *ele, double theta, double baz,
        double lat, double lon, double dep, bool dumpStrain, bool dumpCurl);
    
    // before time loop
    void initialize();
    
    // after time loop
    void finalize();
    
    // record at a time step
    void record(int tstep, double t);
    
    // dump to user-specified format
    void dumpToFile();
    
    // add IO
    void addIO(PointwiseIO *io) {mIOs.push_back(io);};
    
private:
    std::vector<PointwiseInfo> mPointwiseInfo;
    
    // interval
    int mTotalRecordSteps;
    int mRecordInterval;
    
    // buffer
    int mBufferSize;
    int mBufferLine;
    // need row-major to be consistent with netcdf
    RMatXX_RM mBufferDisp;
    RMatXX_RM mBufferStrain;
    RMatXX_RM mBufferCurl;
    RDColX mBufferTime;
    
    // components
    std::string mComponents;
    
    // IO
    std::vector<PointwiseIO *> mIOs;    
    
    // source location
    double mSrcLat, mSrcLon, mSrcDep;
};


