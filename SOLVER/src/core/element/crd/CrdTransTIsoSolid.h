// CrdTransTIsoSolid.h
// created by Kuangdai on 19-May-2017 
// coordinate transformation between (s, phi, z) and (theta, phi, r)

#pragma once

#include "eigenc.h"
#include "eigenp.h"

class CrdTransTIsoSolid {
public:
    
    CrdTransTIsoSolid(const RDMatPP &theta);
    ~CrdTransTIsoSolid() {};
    
    void transformSPZ_RTZ(vec_ar6_CMatPP &u, int Nu) const;
    void transformRTZ_SPZ(vec_ar6_CMatPP &u, int Nu) const;
    
    void transformSPZ_RTZ(vec_ar9_CMatPP &u, int Nu) const;
    void transformRTZ_SPZ(vec_ar9_CMatPP &u, int Nu) const;
    
private:
    
    // t: theta of GLL-points
    RMatPP mSin1t;
    RMatPP mCos1t;
    RMatPP mSin2t;
    RMatPP mCos2t;
};

