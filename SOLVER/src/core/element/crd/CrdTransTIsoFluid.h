// CrdTransTIsoFluid.h
// created by Kuangdai on 19-May-2017 
// coordinate transformation between (s, phi, z) and (theta, phi, r)

#pragma once

#include "eigenc.h"
#include "eigenp.h"

class CrdTransTIsoFluid {
public:
    
    CrdTransTIsoFluid(const RDMatPP &theta);
    ~CrdTransTIsoFluid() {};
    
    void transformSPZ_RTZ(vec_ar3_CMatPP &u, int Nu) const;
    void transformRTZ_SPZ(vec_ar3_CMatPP &u, int Nu) const;
    
private:
    
    // t: theta of GLL-points
    RMatPP mSin1t;
    RMatPP mCos1t;
};

