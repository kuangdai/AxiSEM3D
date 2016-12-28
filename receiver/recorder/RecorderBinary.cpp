// RecorderBinaryAscii.h
// created by Kuangdai on 7-Apr-2016 
// ascii seismogram output

#include "RecorderBinary.h"

RecorderBinary::RecorderBinary(int bufSize, const std::string &fname, bool append):
mBufferSize(bufSize), mFileName(fname), mAppend(append) {
    if (mBufferSize <= 0) mBufferSize = 1;
    mBufferSize = mBufferSize * 4;
    mBufferData.resize(mBufferSize);
    mBufferLoc = 0;
}

void RecorderBinary::open() {
    std::fstream::openmode mode = std::fstream::out | std::fstream::binary;
    if (mAppend) mode = mode | std::fstream::app;
    mFStream.open(mFileName, mode);
    if (!mFStream) throw std::runtime_error("RecorderBinary::open || "
        "Error opening output file: ||" + mFileName);
}

void RecorderBinary::close() {
    if (mFStream.is_open())
        mFStream.close();
}

void RecorderBinary::record(Real t, const RRow3 &u) {
    mBufferData[mBufferLoc++] = t;
    mBufferData[mBufferLoc++] = u(0);
    mBufferData[mBufferLoc++] = u(1);
    mBufferData[mBufferLoc++] = u(2);
    // dump and clear buffer
    if (mBufferLoc == mBufferSize) {
        dumpBufferToFile();
    }
}

void RecorderBinary::dumpBufferToFile() {
    mFStream.write((char *) &(mBufferData[0]), sizeof(Real) * (mBufferLoc));
    mFStream.flush();
    mBufferLoc = 0;
}
