// Station.cpp
// created by Kuangdai on 7-Apr-2016 
// A Station object contains a Seismometer object to compute displacements and 
// a Recorder object to dump displacements.  

#include "Station.h"

Station::Station(int interval, Seismometer *seismometer, Recorder *recorder):
mInterval(interval), mSeismometer(seismometer), mRecorder(recorder) {
    if (mInterval <= 0) mInterval = 1;
    mRecorder->open();
}

Station::~Station() {
    mRecorder->close();
    delete mRecorder;
    delete mSeismometer;
}

void Station::record(int tstep, Real t) {
    if (tstep % mInterval != 0) return;
    static RRow3 gm;
    mSeismometer->getGroundMotion(gm);
    mRecorder->record(t, gm);
}

void Station::dumpLeft() {
    mRecorder->dumpBufferToFile();
}