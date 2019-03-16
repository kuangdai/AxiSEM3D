// SurfaceIO.h
// created by Kuangdai on 28-Nov-2017 
// NetCDF IO for surface wavefield

#pragma once

class NetCDF_Writer;
class SurfaceInfo;
#include "eigenc.h"
#include "eigenp.h"

class SurfaceIO {
public:
    SurfaceIO(bool assemble): mAssemble(assemble) {
        // nothing
    }
    
    // before time loop
    void initialize(int totalRecordSteps, int bufferSize,
        const std::vector<SurfaceInfo> &surfaceInfo,
        double srcLat, double srcLon, double srcDep);
    
    // after time loop
    void finalize();
    
    // dump to netcdf
    void dumpToFile(const std::vector<CMatXX_RM> &bufferDisp, 
        const RDColX &bufferTime, int bufferLine);
    
private:
    // variable names
    std::vector<std::string> mVarNames;
    std::vector<int> mNu;
    
    // file ID
    NetCDF_Writer *mNetCDF = 0;
    
    // location in nc 
    int mCurrentRow = 0;
    
    // minimum MPI rank that has elements
    int mMinRankWithEle = -1;
    
    // source location
    double mSrcLat, mSrcLon, mSrcDep;
    
    // assemble or not
    bool mAssemble = true;
};

