// PointwiseIONetCDF.h
// created by Kuangdai on 1-Jun-2017 
// NetCDF IO for point-wise receivers

#pragma once

#include "PointwiseIO.h"
class NetCDF_Writer;

class PointwiseIONetCDF: public PointwiseIO {
public:
    PointwiseIONetCDF(bool assemble): mAssemble(assemble) {
        // nothing
    }
    
    // before time loop
    void initialize(int totalRecordSteps, int bufferSize, 
        const std::string &components, const std::vector<PointwiseInfo> &receivers, 
        double srcLat, double srcLon, double srcDep);
    
    // after time loop
    void finalize();
    
    // dump to user-specified format
    void dumpToFile(const RMatXX_RM &bufferDisp, 
        const RMatXX_RM &bufferStrain, 
        const RMatXX_RM &bufferCurl, 
        const RDColX &bufferTime, int bufferLine);
    
private:
    // receivers
    const std::vector<PointwiseInfo> *mReceivers;
    
    // variable names
    std::vector<std::string> mVarNamesDisp;
    std::vector<std::string> mVarNamesStrain;
    std::vector<int> mStrainIndex;
    std::vector<std::string> mVarNamesCurl;
    std::vector<int> mCurlIndex;
    
    // file ID
    std::vector<NetCDF_Writer *> mNetCDFs;
    
    // location in nc 
    int mCurrentRow = 0;
    
    // minimum MPI rank that has receivers
    int mMinRankWithRec = -1;
    
    // source location
    double mSrcLat, mSrcLon, mSrcDep;
    
    // assemble or not
    bool mAssemble = true;
    
    // maximum number of stations per file
    const int mMaxNumRecPerFile = 10000;
};

