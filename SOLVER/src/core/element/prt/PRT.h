// PRT.h
// created by Kuangdai on 19-May-2017 
// particle relabelling transformation

#pragma once

class SolidElementResponse;
class FluidElementResponse;
#include <string>
#include "eigenp.h"

class PRT {
public:
    
    virtual ~PRT() {};
    
    virtual void sphericalToUndulated(FluidElementResponse &response) const = 0;
    virtual void undulatedToSpherical(FluidElementResponse &response) const = 0;
    
    virtual void sphericalToUndulated(SolidElementResponse &response) const = 0;
    virtual void undulatedToSpherical(SolidElementResponse &response) const = 0;
    
    virtual std::string verbose() const = 0;
    virtual bool is1D() const = 0;
    
    virtual void checkCompatibility(int Nr) const {};
    
    static PRT *createPRT(const RDMatXN4 &X, bool elem1D);
};
