// Acoustic3X.h
// created by Kuangdai on 6-Jun-2016 
// 3D acoustic with topography

#pragma once

#include "Acoustic.h"
#include "eigenp.h"

class Acoustic3X: public Acoustic {
public:
    // constructor
    Acoustic3X(const RDRowN &theta, const RMatXN KJ, const RMatXN4 X);
    
    // STEP 2: strain ==> stress
    void strainToStress(const vec_ar3_CMatPP &strain, vec_ar3_CMatPP &stress, int Nu) const;
    
    // check size
    void checkCompatibility(int Nr) const;
    
    // verbose
    std::string verbose() const {return "Acoustic3X";};
    
    // change data structure
    // make flat
    static void flattenScalar(const vec_ar3_CMatPP &mat, CMatXN3 &row, int Nu);
    // make structured
    static void stackupScalar(const CMatXN3 &row, vec_ar3_CMatPP &mat, int Nu);

private:

    void rotateStrainToTIso(const CMatXN3 &strainCyln, CMatXN3 &strainTIso, int Nu) const;
    void rotateStressToCyln(const CMatXN3 &stressTIso, CMatXN3 &stressCyln, int Nu) const;
    
    RMatXN mSint, mCost;
    RMatXN mX00, mX01, mX02, mX123;
};

