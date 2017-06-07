// PointwiseIOAscii.cpp
// created by Kuangdai on 1-Jun-2017 
// ascii IO for point-wise receivers

#include "PointwiseIOAscii.h"
#include "Parameters.h"

void PointwiseIOAscii::initialize(int totalRecordSteps, int bufferSize, bool ENZ,
    const std::vector<std::string> &names,
    const std::vector<std::string> &networks) {
    // number
    int numRec = names.size();
    mFileNames.resize(numRec);    
    mFiles.resize(numRec);
    mBuffer = RMatXX::Zero(bufferSize, 4);
    
    // files
    std::string outdir = Parameters::sOutputDirectory + "/stations/";
    for (int irec = 0; irec < numRec; irec++) {
        mFileNames[irec] = outdir + networks[irec] + "_" + names[irec];
        mFileNames[irec] += ENZ ? ".ENZ.ascii" : ".RTZ.ascii";
        mFiles[irec].open(mFileNames[irec], std::fstream::out);
        if (!mFiles[irec]) {
            throw std::runtime_error("PointwiseIOAscii::initialize || "
                "Error opening ascii output file: || " + mFileNames[irec]
                + " || Use NetCDF or ASDF instead of ascii if there are too many stations.");
        }
    }
}

void PointwiseIOAscii::finalize() {
    int numRec = mFileNames.size();
    for (int irec = 0; irec < numRec; irec++) {
        mFiles[irec].close();
    }
}

void PointwiseIOAscii::dumpToFile(const RMatXX &bufferDisp, const RColX &bufferTime, int bufferLine) {
    if (bufferLine == 0) {
        return;
    }
    
    #ifndef NDEBUG
        Eigen::internal::set_is_malloc_allowed(true);
    #endif
    
    int numRec = mFileNames.size();
    for (int irec = 0; irec < numRec; irec++) {
        mBuffer.topRows(bufferLine) << bufferTime.topRows(bufferLine), 
                                       bufferDisp.block(0, irec * 3, bufferLine, 3);
        mFiles[irec] << mBuffer.topRows(bufferLine).format(EIGEN_FMT) << std::endl;   
        mFiles[irec].flush();
    }
    
    #ifndef NDEBUG
        Eigen::internal::set_is_malloc_allowed(false);
    #endif
}

