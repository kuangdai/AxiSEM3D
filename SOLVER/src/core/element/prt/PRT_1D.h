// PRT_1D.h
// created by Kuangdai on 19-May-2017 
// 1D particle relabelling transformation 

#pragma once

#include "PRT.h"
#include "eigenc.h"

class PRT_1D: public PRT {
public:
    PRT_1D(const std::array<RMatPP, 4> &X): mXStruct(X) {};
    ~PRT_1D() {};
    
    void sphericalToUndulated(FluidResponse &response) const;
    void undulatedToSpherical(FluidResponse &response) const;
    
    void sphericalToUndulated(SolidResponse &response) const;
    void undulatedToSpherical(SolidResponse &response) const;
    
    // 9 to 9, for curl computation
    void sphericalToUndulated9(SolidResponse &response) const;
    
    std::string verbose() const {return "PRT_1D";};
    bool is1D() const {return true;};

private:
    std::array<RMatPP, 4> mXStruct;
};

