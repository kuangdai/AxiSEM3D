// Station.h
// created by Kuangdai on 7-Apr-2016 
// A Station object contains a Seismometer object to compute displacements and 
// a Recorder object to dump displacements.  

#pragma once

#include "Seismometer.h"
#include "Recorder.h"

class Station {
public:        
    Station(int interval, Seismometer *seismometer, Recorder *recorder);
    virtual ~Station();
    
    void record(int tstep, Real t);
    
    void dumpLeft();
    
protected:
    int mInterval;
    Seismometer *mSeismometer;
    Recorder *mRecorder;
};