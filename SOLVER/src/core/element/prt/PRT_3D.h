// PRT_3D.h
// created by Kuangdai on 19-May-2017 
// 3D particle relabelling transformation

#pragma once

#include "PRT.h"
#include "eigenc.h"

class PRT_3D: public PRT {
public:
    PRT_3D(const RMatXN4 &X);
    ~PRT_3D() {};
    
    void sphericalToUndulated(FluidResponse &response) const;
    void undulatedToSpherical(FluidResponse &response) const;
    
    void sphericalToUndulated(SolidResponse &response) const;
    void undulatedToSpherical(SolidResponse &response) const;
    
    // 9 to 9, for curl computation
    void sphericalToUndulated9(SolidResponse &response) const;
    
    std::string verbose() const {return "PRT_3D";};
    bool is1D() const {return false;};
    
    void checkCompatibility(int Nr) const;
    
private:
    RMatXN mXFlat0;
    RMatXN mXFlat1;
    RMatXN mXFlat2;
    RMatXN mXFlat3;
};

