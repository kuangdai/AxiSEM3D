// PointwiseRecorder.h
// created by Kuangdai on 1-Jun-2017 
// recorder for point-wise receivers

#pragma once

#include "eigenc.h"
#include "eigenp.h"
class Element;
class PointwiseIO;

class PointwiseRecorder {
public:
    PointwiseRecorder(int totalStepsSTF, int recordInterval, 
        int bufferSize, bool ENZ);
    ~PointwiseRecorder();
    
    // add a receiver
    void addReceiver(const std::string &name, const std::string &network,
        double phi, const RDMatPP &weights, const Element *ele, double theta, double baz);
    
    // before time loop
    void initialize();
    
    // after time loop
    void finalize();
    
    // record at a time step
    void record(int tstep, Real t);
    
    // dump to user-specified format
    void dumpToFile();
    
    // add IO
    void addIO(PointwiseIO *io) {mIOs.push_back(io);};
    
private:
    // receiver info
    struct PointwiseInfo {
        PointwiseInfo(const std::string &name, const std::string &network, 
            double phi, const RDMatPP &weights, const Element *ele,
            double theta, double baz):
            mName(name), mNetwork(network), 
            mPhi(phi), mWeights(weights.cast<Real>()), mElement(ele),
            mTheta(theta), mBAz(baz) {};
        
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
    };
    std::vector<PointwiseInfo> mPointwiseInfo;
    
    // interval
    int mTotalRecordSteps;
    int mRecordInterval;
    
    // buffer
    int mBufferSize;
    int mBufferLine;
    RMatXX mBufferDisp;
    RColX mBufferTime;
    
    // components
    bool mENZ;
    
    // IO
    std::vector<PointwiseIO *> mIOs;    
};


