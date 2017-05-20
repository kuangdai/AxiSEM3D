// PRT_1D.h
// created by Kuangdai on 19-May-2017 
// 1D particle relabelling transformation 

#pragma once

#include "PRT.h"
#include "eigenc.h"

class PRT_1D: public PRT {
public:
    PRT_1D(const RDMatXN4 &X);
    ~PRT_1D() {};
    
    void sphericalToUndulated(FluidElementResponse &response) const;
    void undulatedToSpherical(FluidElementResponse &response) const;
    
    void sphericalToUndulated(SolidElementResponse &response) const;
    void undulatedToSpherical(SolidElementResponse &response) const;
    
    std::string verbose() const {return "PRT_1D";};
    bool is1D() const {return true;};

private:
    std::array<RMatPP, 4> mXStruct;
};

