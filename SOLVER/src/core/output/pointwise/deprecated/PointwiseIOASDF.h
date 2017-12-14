// PointwiseIOASDF.h
// created by Kuangdai on 1-Jun-2017 
// ASDF IO for point-wise receivers

#pragma once

#include "PointwiseIO.h"
#include "boost/date_time/posix_time/posix_time_types.hpp"
class NetCDF_Writer;

class PointwiseIOASDF: public PointwiseIO {
public:
    PointwiseIOASDF(double srcLat, double srcLon, double srcDep,
        const std::string &sourceFile):
        mSrcLat(srcLat), mSrcLon(srcLon), mSrcDep(srcDep), 
        mSourceFile(sourceFile) {};
    
    // before time loop
    void initialize(int totalRecordSteps, int bufferSize, 
        const std::string &components, const std::vector<PointwiseInfo> &receivers);
    
    // after time loop
    void finalize();
    
    // dump to user-specified format
    void dumpToFile(const RMatXX_RM &bufferDisp, const RColX &bufferTime, int bufferLine);
    
private:
    void createQuakeML(NetCDF_Writer &nw, std::string &sourceID, std::string &sourceT0_UTC);
    void createStationML(NetCDF_Writer &nw, int irec, const std::string &sourceID);
    
    static boost::posix_time::ptime UTCfromString(const std::string &utcStr);
    static std::string UTCToString(const boost::posix_time::ptime &utc, bool printfrac);
    
    // variable names
    std::vector<std::string> mVarNames;
    
    // file ID
    NetCDF_Writer *mNetCDF = 0;
    
    // location in nc 
    int mCurrentRow = 0;
    
    // header info
    std::string mComponents;
    double mSrcLat;
    double mSrcLon;
    double mSrcDep;
    std::string mSourceFile;
    std::vector<std::string> mNetworks;
    std::vector<std::string> mNames;
    std::vector<double> mLats;
    std::vector<double> mLons;
    std::vector<double> mDeps;
};

