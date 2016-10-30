// RecorderAsciiAscii.h
// created by Kuangdai on 7-Apr-2016 
// ascii seismogram output

#include "RecorderAscii.h"

RecorderAscii::RecorderAscii(int bufSize, const std::string &fname, bool append):
mBufferSize(bufSize), mBufferLine(0), mFileName(fname), mAppend(append) {
    if (mBufferSize <= 0) mBufferSize = 1;
    mBufferData = Eigen::Matrix<Real, Eigen::Dynamic, 4>(mBufferSize, 4);
}

void RecorderAscii::open() {
    std::fstream::openmode mode = std::fstream::out;
    if (mAppend) mode = mode | std::fstream::app;
    mFStream.open(mFileName, mode);
    if (!mFStream) throw std::runtime_error("RecorderAscii::open || "
        "Error opening output file: ||" + mFileName);
}

void RecorderAscii::close() {
    if (mFStream.is_open())
        mFStream.close();
}

void RecorderAscii::record(Real t, const RRow3 &u) {
    mBufferData(mBufferLine, 0) = t;
    mBufferData.block(mBufferLine, 1, 1, 3) = u;
    mBufferLine++;
    // dump and clear buffer
    if (mBufferLine == mBufferSize) {
        dumpBufferToFile();
    }
}

void RecorderAscii::dumpBufferToFile() {
    // each row starts with a space
    #ifndef NDEBUG
        Eigen::internal::set_is_malloc_allowed(true);
        mFStream << mBufferData.topRows(mBufferLine).format(EIGEN_FMT) << std::endl;   
        Eigen::internal::set_is_malloc_allowed(false); 
    #else
        mFStream << mBufferData.topRows(mBufferLine).format(EIGEN_FMT) << std::endl;  
    #endif
    mFStream.flush();
    mBufferLine = 0;
}
