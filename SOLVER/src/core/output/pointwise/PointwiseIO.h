// PointwiseIO.h
// created by Kuangdai on 1-Jun-2017 
// IO for point-wise receivers

#pragma once

#include "eigenc.h"
#include "eigenp.h"
#include <string>
class PointwiseInfo;

class PointwiseIO {
public:
    virtual ~PointwiseIO() {};
    
    // before time loop
    virtual void initialize(int totalRecordSteps, int bufferSize, 
        const std::string &components, const std::vector<PointwiseInfo> &receivers,
        double srcLat, double srcLon, double srcDep) = 0;
    
    // after time loop
    virtual void finalize() = 0;
    
    // dump to user-specified format
    virtual void dumpToFile(const RMatXX_RM &bufferDisp, 
        const RMatXX_RM &bufferStrain, 
        const RMatXX_RM &bufferCurl, 
        const RDColX &bufferTime, int bufferLine) = 0;
};

