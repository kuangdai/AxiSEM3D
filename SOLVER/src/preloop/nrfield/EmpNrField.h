// EmpNrField.h
// created by Kuangdai on 13-May-2016 
// empirical nr integer field

#pragma once
#include "NrField.h"

class EmpNrField: public NrField {
public:
    EmpNrField(bool useLucky, int nu_ref, int nu_min, 
        bool scaleS = true, bool scaleT = true, bool scaleD = true, 
        double powS = 1.,
        double factPI = 5., double startT = pi / 4, double powT = 3.,
        double factD0 = 2., double startD = 200e3, double endD = 300e3);
    
    int getNrAtPoint(const RDCol2 &coords) const;
    
    std::string verbose() const;
    
private:
    // basic values
    int mNuRef, mNuMin;
    
    // distance to axis
    // (s / mROuter) ^ mPowS
    bool mScaleS;
    double mROuter, mPowS;
    
    // epicentral distance
    // 1. + (mFactPI - 1.) * ((theta - mStartT) / (pi - mStartT)) ^ mPowT
    bool mScaleT;
    double mFactPI, mStartT, mPowT;
    
    // surface wave 
    // d = 0 (surface)  -> factor = mFactD0
    // d = mStartD      -> factor = mFactD0
    // d = mEndD        -> factor = 1.0 
    bool mScaleD;
    double mFactD0, mStartD, mEndD; 
};

