// STF.h
// created by Kuangdai on 9-May-2016
// source time function

#pragma once

#include <vector>
#include <string>
#include "global.h"

class Domain;
class Parameters;

class STF {
public:
    virtual ~STF() {};

    void release(Domain &domain) const;

    virtual std::string verbose() const = 0;

    static void buildInparam(STF *&stf, const Parameters &par, double dt, int verbose);
    
    int getSize() const {return mSTF.size();};

protected:
    double mDeltaT;
    double mShift;
    std::vector<double> mSTF;
};
