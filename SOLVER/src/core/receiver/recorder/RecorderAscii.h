// RecorderAscii.h
// created by Kuangdai on 7-Apr-2016 
// ascii seismogram output

#pragma once

#include "Recorder.h"
#include <fstream>

class RecorderAscii: public Recorder {
public: 
    RecorderAscii(int bufSize, const std::string &fname, bool append);
    
    void open();
    void close();
    void record(Real t, const RRow3 &u);
    void dumpBufferToFile();

protected:
    // buffer size
    int mBufferSize;
    // current line in buffer
    int mBufferLine;
    // data buffer
    Eigen::Matrix<Real, Eigen::Dynamic, 4, Eigen::RowMajor> mBufferData;
    
    // file name
    std::string mFileName;
    // file options
    bool mAppend;
    // file stream
    std::fstream mFStream;
    
};
