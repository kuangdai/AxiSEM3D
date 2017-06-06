// Receiver.h
// created by Kuangdai on 1-Jun-2016 
// receiver collections

#pragma once

#include <vector>
#include <string>

class Domain;
class Mesh;
class Parameters;
class Receiver;

class ReceiverCollection {
public:
    ReceiverCollection(const std::string &fileRec, bool geographic, 
        double srcLat, double srcLon, double srcDep);
    ~ReceiverCollection();
    
    void release(Domain &domain, const Mesh &mesh); 
    
    std::string verbose() const;
    
    static void buildInparam(ReceiverCollection *&rec, const Parameters &par, 
        double srcLat, double srcLon, double srcDep, int totalStepsSTF, int verbose);
        
private:
    
    // receivers
    std::vector<Receiver *> mReceivers;
    
    // input
    std::string mInputFile;
    bool mGeographic;
    
    // options
    int mTotalStepsSTF = 0;
    int mRecordInterval = 1;
    int mBufferSize = 1000;
    bool mENZ = false;
    
    // IO
    bool mAscii = false;
    bool mNetCDF = false;
    bool mASDF = false;
    
    // for verbose
    int mWidthName;
    int mWidthNetwork;
};

