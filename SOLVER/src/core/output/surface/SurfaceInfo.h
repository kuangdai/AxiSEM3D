// SurfaceInfo.h
// created by Kuangdai on 27-Nov-2017 
// surface element info

#pragma once

class Element;
#include "eigenc.h"

class SurfaceInfo {
public:
    SurfaceInfo(const Element *ele, int surfSide);
    ~SurfaceInfo() {};
    
    void setGlobalTag(int globalTag) {mGlobalTag = globalTag;};
    int getGlobalTag() const {return mGlobalTag;};
    
    double getThetaMin() const {return std::min(mTheta0, mTheta1);};
    double getTheta0() const {return mTheta0;};
    double getTheta1() const {return mTheta1;};
    
    void initBuffer(int bufferSize, CMatXX_RM &bufferDisp);
    void feedBuffer(int bufferLine, CMatXX_RM &bufferDisp);
    
    int getMaxNu() const;
    
private:    
    const Element *mElement;
    int mSurfSide;
    int mGlobalTag = -1;
    
    // redundant but useful
    double mTheta0;
    double mTheta1;
};

