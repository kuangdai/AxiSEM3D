// PointwiseIOASDF.h
// created by Kuangdai on 1-Jun-2017 
// ASDF IO for point-wise receivers

#pragma once

#include "PointwiseIO.h"
class NetCDF_Writer;

class PointwiseIOASDF: public PointwiseIO {
public:
    PointwiseIOASDF(double srcLat, double srcLon, double srcDep):
        mSrcLat(srcLat), mSrcLon(srcLon), mSrcDep(srcDep) {};
    
    // before time loop
    void initialize(int totalRecordSteps, int bufferSize, bool ENZ, 
        const std::vector<PointwiseInfo> &receivers);
    
    // after time loop
    void finalize();
    
    // dump to user-specified format
    void dumpToFile(const RMatXX_RM &bufferDisp, const RColX &bufferTime, int bufferLine);
    
private:
    void createQuakeML(NetCDF_Writer &nw);
    void createStationML(NetCDF_Writer &nw, int irec);
    
    
    // variable names
    std::vector<std::string> mVarNames;
    
    // file ID
    NetCDF_Writer *mNetCDF = 0;
    
    // location in nc 
    int mCurrentRow = 0;
    
    // header info
    bool mENZ;
    double mSrcLat;
    double mSrcLon;
    double mSrcDep;
    std::vector<std::string> mNetworks;
    std::vector<std::string> mNames;
    std::vector<double> mLats;
    std::vector<double> mLons;
    std::vector<double> mDeps;
};

