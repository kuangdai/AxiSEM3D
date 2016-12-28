// Recorder.h
// created by Kuangdai on 7-Apr-2016 
// seismogram output
// implement a new sub-class for other formats

#pragma once

#include "global.h"
#include "eigenc.h"

class Recorder {
public: 
    virtual ~Recorder() {};
    virtual void open() = 0;
    virtual void close() = 0;
    virtual void record(Real t, const RRow3 &u) = 0;
    virtual void dumpBufferToFile() = 0;
};
