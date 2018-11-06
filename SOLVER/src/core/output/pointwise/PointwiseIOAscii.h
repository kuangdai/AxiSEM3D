// PointwiseIOAscii.h
// created by Kuangdai on 1-Jun-2017 
// ascii IO for point-wise receivers

#pragma once

#include <fstream>
#include "PointwiseIO.h"

class PointwiseIOAscii: public PointwiseIO {
public:
    // before time loop
    void initialize(int totalRecordSteps, int bufferSize, 
        const std::string &components, const std::vector<PointwiseInfo> &receivers, 
        double srcLat, double srcLon, double srcDep);
    
    // after time loop
    void finalize();
    
    // dump to user-specified format
    void dumpToFile(const RMatXX_RM &bufferDisp, 
        const RMatXX_RM &bufferStrain,
        const RMatXX_RM &bufferCurl,
        const RDColX &bufferTime, int bufferLine);
    
private:
    // file names
    std::vector<std::string> mFileNamesDisp;
    
    // fstream
    std::vector<std::fstream *> mFilesDisp;
    
    // buffer
    RDMatXX mBufferDisp;
    
    // file names
    std::vector<std::string> mFileNamesStrain;
    
    // fstream
    std::vector<std::fstream *> mFilesStrain;
    
    // buffer
    RDMatXX mBufferStrain;
    
    // file names
    std::vector<std::string> mFileNamesCurl;
    
    // fstream
    std::vector<std::fstream *> mFilesCurl;
    
    // buffer
    RDMatXX mBufferCurl;
};

