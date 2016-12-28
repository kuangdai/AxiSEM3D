// RecorderBinary.h
// created by Kuangdai on 7-Apr-2016 
// binary seismogram output

#pragma once

#include "Recorder.h"
#include <fstream>
#include <vector>

class RecorderBinary: public Recorder
{
public: 
    RecorderBinary(int bufSize, const std::string &fname, bool append);
    
    void open();
    void close();
    void record(Real t, const RRow3 &u);
    void dumpBufferToFile();

protected:
    // total buffer size
    int mBufferSize;
    // current position in buffer
    int mBufferLoc;
    // data buffer
    std::vector<Real> mBufferData;
    
    // file name
    std::string mFileName;
    // file options
    bool mAppend;
    // file stream
    std::fstream mFStream;
    
};
