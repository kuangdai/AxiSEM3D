// Mass.h
// created by Kuangdai on 3-Apr-2016 
// base class of mass

#pragma once

#include "eigenc.h"

class Mass {
public:
    virtual ~Mass() {};
    
    // compute accel in-place
    virtual void computeAccel(CMatX3 &stiff) const = 0;
    virtual void computeAccel(CColX &stiff) const = 0;
    
    // check compatibility
    virtual void checkCompatibility(int nr) const {};
    
    // verbose 
    virtual std::string verbose() const = 0;
};
