// Elastic3D.h
// created by Kuangdai on 29-Apr-2016 
// base class of 3D elasticity

#pragma once

#include "eigenc.h"
#include "Elastic.h"

class Attenuation3D;

class Elastic3D: public Elastic {
public:
    Elastic3D(Attenuation3D *att);
    virtual ~Elastic3D();
    
    // check compatibility
    virtual void checkCompatibility(int Nr, bool isVoigt) const; 
    
    // reset to zero 
    void resetZero(); 
    
    // change data structure
    // make flat
    static void flattenVector(const vec_ar9_CMatPP &mat, CMatXN9 &row, int Nu);
    static void flattenVectorVoigt(const vec_ar9_CMatPP &mat, CMatXN6 &row, int Nu);
    // make structured
    static void stackupVector(const CMatXN9 &row, vec_ar9_CMatPP &mat, int Nu);
    static void stackupVectorVoigt(const CMatXN6 &row, vec_ar9_CMatPP &mat, int Nu);
    
protected:
    Attenuation3D *mAttenuation;
};
