// PointwiseIOAscii.h
// created by Kuangdai on 1-Jun-2017 
// ascii IO for point-wise receivers

#pragma once

#include <fstream>
#include "PointwiseIO.h"

class PointwiseIOAscii: public PointwiseIO {
public:
    // before time loop
    void initialize(int totalRecordSteps, int bufferSize, bool ENZ, 
        const std::vector<PointwiseInfo> &receivers);
    
    // after time loop
    void finalize();
    
    // dump to user-specified format
    void dumpToFile(const RMatXX_RM &bufferDisp, const RColX &bufferTime, int bufferLine);
    
private:
    // file names
    std::vector<std::string> mFileNames;
    
    // fstream
    std::vector<std::fstream> mFiles;
    
    // buffer
    RMatXX mBuffer;
};

