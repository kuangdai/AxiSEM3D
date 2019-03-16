// SurfaceRecorder.h
// created by Kuangdai on 27-Nov-2017 
// recorder for surface wavefield

#pragma once

#include "eigenc.h"
#include "eigenp.h"
class SurfaceInfo;
class SurfaceIO;
class Element;

class SurfaceRecorder {
public:
    SurfaceRecorder(int totalRecordSteps, int recordInterval, int bufferSize,
        double srcLat, double srcLon, double srcDep, bool assemble);
    ~SurfaceRecorder();

    // add a surface element
    void addElement(Element *ele, int surfSide);

    // before time loop
    void initialize();

    // after time loop
    void finalize();

    // record at a time step
    void record(int tstep, double t);

    // dump to netcdf
    void dumpToFile();

private:
    // surface elements
    std::vector<SurfaceInfo> mSurfaceInfo;
    
    // interval
    int mTotalRecordSteps;
    int mRecordInterval;

    // buffer
    int mBufferSize;
    int mBufferLine;
    
    // buffer
    RDColX mBufferTime;
    std::vector<CMatXX_RM> mBufferDisp;
    
    // IO
    SurfaceIO *mIO;
    
    // source location
    double mSrcLat, mSrcLon, mSrcDep;
};


